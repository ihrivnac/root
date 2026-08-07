// Author: ROOT team

#include "TGeoMultiUnion.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TRandom3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <numeric>
#include <typeinfo>

ClassImp(TGeoMultiUnion);

namespace {
constexpr Int_t kBoxStride = 6;
constexpr Int_t kSmallUnionNodeLimit = 2;
constexpr std::size_t kBVHStackSize = 64;
constexpr Int_t kNodeParameterStride = 6;
constexpr UChar_t kPrimitiveMask = 0x3;
constexpr UChar_t kGenericNode = 0;
constexpr UChar_t kBoxNode = 1;
constexpr UChar_t kTubeNode = 2;
constexpr UChar_t kRotatedNode = 0x4;

Bool_t ContainsBox(const Double_t *box, const Double_t *point)
{
   const Double_t tolerance = TGeoShape::Tolerance();
   return point[0] >= box[0] - tolerance && point[0] <= box[1] + tolerance &&
          point[1] >= box[2] - tolerance && point[1] <= box[3] + tolerance &&
          point[2] >= box[4] - tolerance && point[2] <= box[5] + tolerance;
}

Double_t BoxSafety(const Double_t *box, const Double_t *point)
{
   Double_t result = 0.;
   const Double_t tolerance = TGeoShape::Tolerance();
   for (Int_t axis = 0; axis < 3; ++axis) {
      const Double_t gap = std::max(box[2 * axis] - point[axis], point[axis] - box[2 * axis + 1]);
      result = std::max(result, gap - tolerance);
   }
   return result;
}

Double_t BoxRayDistance(const Double_t *box, const Double_t *point, const Double_t *dir, Double_t step)
{
   Double_t enter = 0.;
   Double_t leave = step;
   for (Int_t axis = 0; axis < 3; ++axis) {
      const Double_t lo = box[2 * axis];
      const Double_t hi = box[2 * axis + 1];
      if (std::abs(dir[axis]) < std::numeric_limits<Double_t>::epsilon()) {
         if (point[axis] < lo || point[axis] > hi)
            return TGeoShape::Big();
         continue;
      }
      Double_t first = (lo - point[axis]) / dir[axis];
      Double_t last = (hi - point[axis]) / dir[axis];
      if (first > last)
         std::swap(first, last);
      enter = std::max(enter, first);
      leave = std::min(leave, last);
      if (enter > leave)
         return TGeoShape::Big();
   }
   return leave >= 0. ? enter : TGeoShape::Big();
}

Double_t PushDistance(Double_t distance)
{
   return 10. * TGeoShape::Tolerance() * std::max(1., std::abs(distance));
}
} // namespace

TGeoMultiUnion::TGeoMultiUnion() : TGeoBBox(0., 0., 0.)
{
   fMatrices.SetOwner(kTRUE);
}

TGeoMultiUnion::TGeoMultiUnion(const char *name) : TGeoBBox(name, 0., 0., 0.)
{
   fMatrices.SetOwner(kTRUE);
}

TGeoMultiUnion::~TGeoMultiUnion() = default;

void TGeoMultiUnion::ClearThreadData() const
{
   for (Int_t inode = 0; inode < GetNnodes(); ++inode)
      GetShape(inode)->ClearThreadData();
}

void TGeoMultiUnion::CreateThreadData(Int_t nthreads)
{
   for (Int_t inode = 0; inode < GetNnodes(); ++inode)
      GetShape(inode)->CreateThreadData(nthreads);
}

void TGeoMultiUnion::AddNode(TGeoShape &shape, const TGeoMatrix &matrix)
{
   AddNode(&shape, &matrix);
}

void TGeoMultiUnion::AddNode(TGeoShape *shape, const TGeoMatrix *matrix)
{
   if (!shape) {
      Error("AddNode", "Cannot add a null shape to %s", GetName());
      return;
   }
   if (!dynamic_cast<TGeoBBox *>(shape)) {
      Error("AddNode", "Shape %s has no finite bounding box", shape->GetName());
      return;
   }
   if (shape->TestShapeBit(TGeoShape::kGeoHalfSpace)) {
      Error("AddNode", "Infinite half-space %s cannot be a multi-union node", shape->GetName());
      return;
   }
   if (matrix && matrix->IsScale()) {
      Error("AddNode", "Scaled transformations are not supported for node %s", shape->GetName());
      return;
   }

   fShapes.Add(shape);
   fMatrices.Add(matrix ? new TGeoHMatrix(*matrix) : new TGeoHMatrix());
   fVoxelized = kFALSE;
   fHasBVH = kFALSE;
   fBoxes.clear();
   fBVHBoxes.clear();
   fBVHChildren.clear();
   fNodeFlags.clear();
   fNodeTranslations.clear();
   fNodeParameters.clear();
   ComputeBBox();
}

void TGeoMultiUnion::AfterStreamer()
{
   fMatrices.SetOwner(kTRUE);
   if (fVoxelized)
      Voxelize();
   else
      ComputeBBox();
}

TGeoShape *TGeoMultiUnion::GetShape(Int_t inode) const
{
   return inode >= 0 && inode < GetNnodes() ? static_cast<TGeoShape *>(fShapes.UncheckedAt(inode)) : nullptr;
}

TGeoMatrix *TGeoMultiUnion::GetMatrix(Int_t inode) const
{
   return inode >= 0 && inode < GetNnodes() ? static_cast<TGeoMatrix *>(fMatrices.UncheckedAt(inode)) : nullptr;
}

void TGeoMultiUnion::TransformToNode(Int_t inode, const Double_t *point, const Double_t *dir, Double_t *local,
                                     Double_t *localDir) const
{
   TransformPointToNode(inode, point, local);
   if (dir)
      TransformDirectionToNode(inode, dir, localDir);
}

void TGeoMultiUnion::TransformPointToNode(Int_t inode, const Double_t *point, Double_t *local) const
{
   if (fNodeFlags.size() == static_cast<std::size_t>(GetNnodes()) && !(fNodeFlags[inode] & kRotatedNode)) {
      const Double_t *translation = &fNodeTranslations[3 * inode];
      local[0] = point[0] - translation[0];
      local[1] = point[1] - translation[1];
      local[2] = point[2] - translation[2];
      return;
   }
   GetMatrix(inode)->MasterToLocal(point, local);
}

void TGeoMultiUnion::TransformDirectionToNode(Int_t inode, const Double_t *dir, Double_t *localDir) const
{
   if (fNodeFlags.size() == static_cast<std::size_t>(GetNnodes()) && !(fNodeFlags[inode] & kRotatedNode)) {
      std::memcpy(localDir, dir, 3 * sizeof(Double_t));
      return;
   }
   GetMatrix(inode)->MasterToLocalVect(dir, localDir);
}

Bool_t TGeoMultiUnion::NodeContains(Int_t inode, const Double_t *local) const
{
   if (fNodeFlags.size() != static_cast<std::size_t>(GetNnodes()))
      return GetShape(inode)->Contains(local);
   const Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
   switch (fNodeFlags[inode] & kPrimitiveMask) {
   case kBoxNode:
      return std::abs(local[0] - parameters[3]) <= parameters[0] &&
             std::abs(local[1] - parameters[4]) <= parameters[1] &&
             std::abs(local[2] - parameters[5]) <= parameters[2];
   case kTubeNode: {
      if (std::abs(local[2]) > parameters[2])
         return kFALSE;
      const Double_t radiusSquared = local[0] * local[0] + local[1] * local[1];
      return radiusSquared >= parameters[0] * parameters[0] && radiusSquared <= parameters[1] * parameters[1];
   }
   default:
      return GetShape(inode)->Contains(local);
   }
}

Double_t TGeoMultiUnion::NodeSafety(Int_t inode, const Double_t *local, Bool_t in) const
{
   if (fNodeFlags.size() != static_cast<std::size_t>(GetNnodes()))
      return GetShape(inode)->Safety(local, in);
   const Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
   switch (fNodeFlags[inode] & kPrimitiveMask) {
   case kBoxNode: {
      const Double_t sx = parameters[0] - std::abs(local[0] - parameters[3]);
      const Double_t sy = parameters[1] - std::abs(local[1] - parameters[4]);
      const Double_t sz = parameters[2] - std::abs(local[2] - parameters[5]);
      return in ? std::min({sx, sy, sz}) : std::max({-sx, -sy, -sz});
   }
   case kTubeNode:
      return TGeoTube::SafetyS(local, in, parameters[0], parameters[1], parameters[2]);
   default:
      return GetShape(inode)->Safety(local, in);
   }
}

Double_t TGeoMultiUnion::NodeDistFromInside(Int_t inode, const Double_t *local, const Double_t *localDir) const
{
   if (fNodeFlags.size() != static_cast<std::size_t>(GetNnodes()))
      return GetShape(inode)->DistFromInside(local, localDir, 3);
   const Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
   switch (fNodeFlags[inode] & kPrimitiveMask) {
   case kBoxNode:
      return TGeoBBox::DistFromInside(local, localDir, parameters[0], parameters[1], parameters[2], parameters + 3);
   case kTubeNode:
      return TGeoTube::DistFromInsideS(local, localDir, parameters[0], parameters[1], parameters[2]);
   default:
      return GetShape(inode)->DistFromInside(local, localDir, 3);
   }
}

Double_t TGeoMultiUnion::NodeDistFromOutside(Int_t inode, const Double_t *local, const Double_t *localDir,
                                             Double_t step) const
{
   if (fNodeFlags.size() != static_cast<std::size_t>(GetNnodes()))
      return GetShape(inode)->DistFromOutside(local, localDir, 3, step);
   const Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
   switch (fNodeFlags[inode] & kPrimitiveMask) {
   case kBoxNode:
      return TGeoBBox::DistFromOutside(local, localDir, parameters[0], parameters[1], parameters[2], parameters + 3,
                                       step);
   case kTubeNode:
      return TGeoTube::DistFromOutsideS(local, localDir, parameters[0], parameters[1], parameters[2]);
   default:
      return GetShape(inode)->DistFromOutside(local, localDir, 3, step);
   }
}

Bool_t TGeoMultiUnion::AcceptNode(Int_t inode, const Double_t *point) const
{
   if (!fVoxelized || fBoxes.size() != static_cast<std::size_t>(kBoxStride * GetNnodes()))
      return kTRUE;
   const Double_t *box = &fBoxes[kBoxStride * inode];
   const Double_t tol = TGeoShape::Tolerance();
   return point[0] >= box[0] - tol && point[0] <= box[1] + tol && point[1] >= box[2] - tol &&
          point[1] <= box[3] + tol && point[2] >= box[4] - tol && point[2] <= box[5] + tol;
}

void TGeoMultiUnion::BuildBVH()
{
   fHasBVH = kFALSE;
   fBVHBoxes.clear();
   fBVHChildren.clear();
   const Int_t nodeCount = GetNnodes();
   if (nodeCount <= kSmallUnionNodeLimit ||
       fBoxes.size() != static_cast<std::size_t>(kBoxStride * nodeCount))
      return;

   std::vector<Int_t> indices(nodeCount);
   std::iota(indices.begin(), indices.end(), 0);
   fBVHBoxes.reserve(kBoxStride * (2 * nodeCount - 1));
   fBVHChildren.reserve(2 * (2 * nodeCount - 1));

   std::function<Int_t(Int_t, Int_t)> build = [&](Int_t begin, Int_t end) {
      const Int_t treeNode = static_cast<Int_t>(fBVHChildren.size() / 2);
      fBVHChildren.resize(fBVHChildren.size() + 2, -1);
      fBVHBoxes.resize(fBVHBoxes.size() + kBoxStride);
      if (end - begin == 1) {
         fBVHChildren[2 * treeNode + 1] = indices[begin];
         std::copy_n(&fBoxes[kBoxStride * indices[begin]], kBoxStride, &fBVHBoxes[kBoxStride * treeNode]);
         return treeNode;
      }

      Double_t centroidExtent[6] = {TGeoShape::Big(), -TGeoShape::Big(), TGeoShape::Big(),
                                    -TGeoShape::Big(), TGeoShape::Big(), -TGeoShape::Big()};
      for (Int_t position = begin; position < end; ++position) {
         const Double_t *box = &fBoxes[kBoxStride * indices[position]];
         for (Int_t axis = 0; axis < 3; ++axis) {
            const Double_t center = box[2 * axis] + box[2 * axis + 1];
            centroidExtent[2 * axis] = std::min(centroidExtent[2 * axis], center);
            centroidExtent[2 * axis + 1] = std::max(centroidExtent[2 * axis + 1], center);
         }
      }
      Int_t splitAxis = 0;
      for (Int_t axis = 1; axis < 3; ++axis) {
         if (centroidExtent[2 * axis + 1] - centroidExtent[2 * axis] >
             centroidExtent[2 * splitAxis + 1] - centroidExtent[2 * splitAxis])
            splitAxis = axis;
      }
      const Int_t middle = begin + (end - begin) / 2;
      std::nth_element(indices.begin() + begin, indices.begin() + middle, indices.begin() + end,
                       [&](Int_t first, Int_t second) {
                          const Double_t *firstBox = &fBoxes[kBoxStride * first];
                          const Double_t *secondBox = &fBoxes[kBoxStride * second];
                          return firstBox[2 * splitAxis] + firstBox[2 * splitAxis + 1] <
                                 secondBox[2 * splitAxis] + secondBox[2 * splitAxis + 1];
                       });
      const Int_t left = build(begin, middle);
      const Int_t right = build(middle, end);
      fBVHChildren[2 * treeNode] = left;
      fBVHChildren[2 * treeNode + 1] = right;
      Double_t *treeBox = &fBVHBoxes[kBoxStride * treeNode];
      const Double_t *leftBox = &fBVHBoxes[kBoxStride * left];
      const Double_t *rightBox = &fBVHBoxes[kBoxStride * right];
      for (Int_t axis = 0; axis < 3; ++axis) {
         treeBox[2 * axis] = std::min(leftBox[2 * axis], rightBox[2 * axis]);
         treeBox[2 * axis + 1] = std::max(leftBox[2 * axis + 1], rightBox[2 * axis + 1]);
      }
      return treeNode;
   };
   build(0, nodeCount);
   fHasBVH = kTRUE;
}

Bool_t TGeoMultiUnion::CrossesNodeBox(Int_t inode, const Double_t *point, const Double_t *dir, Double_t step) const
{
   if (!fVoxelized || fBoxes.size() != static_cast<std::size_t>(kBoxStride * GetNnodes()))
      return kTRUE;

   return BoxRayDistance(&fBoxes[kBoxStride * inode], point, dir, step) < TGeoShape::Big();
}

void TGeoMultiUnion::Voxelize()
{
   fBoxes.assign(kBoxStride * GetNnodes(), 0.);
   fNodeFlags.assign(GetNnodes(), kGenericNode);
   fNodeTranslations.assign(3 * GetNnodes(), 0.);
   fNodeParameters.assign(kNodeParameterStride * GetNnodes(), 0.);
   Double_t vertices[24];
   Double_t master[3];
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      auto *shape = GetShape(inode);
      const auto *matrix = GetMatrix(inode);
      UChar_t flags = matrix->IsRotation() ? kRotatedNode : 0;
      if (typeid(*shape) == typeid(TGeoBBox)) {
         const auto *box = static_cast<TGeoBBox *>(shape);
         flags |= kBoxNode;
         Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
         parameters[0] = box->GetDX();
         parameters[1] = box->GetDY();
         parameters[2] = box->GetDZ();
         std::copy_n(box->GetOrigin(), 3, parameters + 3);
      } else if (typeid(*shape) == typeid(TGeoTube)) {
         const auto *tube = static_cast<TGeoTube *>(shape);
         flags |= kTubeNode;
         Double_t *parameters = &fNodeParameters[kNodeParameterStride * inode];
         parameters[0] = tube->GetRmin();
         parameters[1] = tube->GetRmax();
         parameters[2] = tube->GetDz();
      }
      fNodeFlags[inode] = flags;
      std::copy_n(matrix->GetTranslation(), 3, &fNodeTranslations[3 * inode]);
      auto *boxShape = static_cast<TGeoBBox *>(shape);
      if (boxShape->IsNullBox())
         shape->ComputeBBox();
      boxShape->SetBoxPoints(vertices);
      Double_t *box = &fBoxes[kBoxStride * inode];
      box[0] = box[2] = box[4] = TGeoShape::Big();
      box[1] = box[3] = box[5] = -TGeoShape::Big();
      for (Int_t vertex = 0; vertex < 8; ++vertex) {
         GetMatrix(inode)->LocalToMaster(&vertices[3 * vertex], master);
         for (Int_t axis = 0; axis < 3; ++axis) {
            box[2 * axis] = std::min(box[2 * axis], master[axis]);
            box[2 * axis + 1] = std::max(box[2 * axis + 1], master[axis]);
         }
      }
   }
   fVoxelized = kTRUE;
   BuildBVH();
   ComputeBBox();
}

void TGeoMultiUnion::ComputeBBox()
{
   if (GetNnodes() == 0) {
      SetBoxDimensions(0., 0., 0.);
      return;
   }

   const Bool_t restoreVoxelized = fVoxelized;
   const Bool_t restoreHasBVH = fHasBVH;
   if (fBoxes.size() != static_cast<std::size_t>(kBoxStride * GetNnodes())) {
      fVoxelized = kFALSE;
      Voxelize();
      fVoxelized = restoreVoxelized;
      fHasBVH = restoreHasBVH;
   }
   Double_t extent[6] = {TGeoShape::Big(), -TGeoShape::Big(), TGeoShape::Big(),
                         -TGeoShape::Big(), TGeoShape::Big(), -TGeoShape::Big()};
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      const Double_t *box = &fBoxes[kBoxStride * inode];
      for (Int_t axis = 0; axis < 3; ++axis) {
         extent[2 * axis] = std::min(extent[2 * axis], box[2 * axis]);
         extent[2 * axis + 1] = std::max(extent[2 * axis + 1], box[2 * axis + 1]);
      }
   }
   Double_t origin[3];
   for (Int_t axis = 0; axis < 3; ++axis)
      origin[axis] = .5 * (extent[2 * axis] + extent[2 * axis + 1]);
   SetBoxDimensions(.5 * (extent[1] - extent[0]), .5 * (extent[3] - extent[2]), .5 * (extent[5] - extent[4]),
                    origin);
}

Bool_t TGeoMultiUnion::Contains(const Double_t *point) const
{
   Double_t local[3];
   if (HasBVH()) {
      std::array<Int_t, kBVHStackSize> stack;
      std::size_t stackSize = 1;
      stack[0] = 0;
      while (stackSize) {
         const Int_t treeNode = stack[--stackSize];
         if (!ContainsBox(&fBVHBoxes[kBoxStride * treeNode], point))
            continue;
         const Int_t left = fBVHChildren[2 * treeNode];
         const Int_t right = fBVHChildren[2 * treeNode + 1];
         if (left < 0) {
            TransformPointToNode(right, point, local);
            if (NodeContains(right, local))
               return kTRUE;
            continue;
         }
         stack[stackSize++] = right;
         stack[stackSize++] = left;
      }
      return kFALSE;
   }
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      if (!AcceptNode(inode, point))
         continue;
      TransformPointToNode(inode, point, local);
      if (NodeContains(inode, local))
         return kTRUE;
   }
   return kFALSE;
}

void TGeoMultiUnion::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
{
   for (Int_t i = 0; i < vecsize; ++i)
      inside[i] = Contains(&points[3 * i]);
}

Double_t TGeoMultiUnion::Safety(const Double_t *point, Bool_t in) const
{
   Double_t result = TGeoShape::Big();
   Double_t local[3];
   if (HasBVH()) {
      std::array<Int_t, kBVHStackSize> stack;
      std::size_t stackSize = 1;
      stack[0] = 0;
      while (stackSize) {
         const Int_t treeNode = stack[--stackSize];
         const Double_t *treeBox = &fBVHBoxes[kBoxStride * treeNode];
         if ((in && !ContainsBox(treeBox, point)) || (!in && BoxSafety(treeBox, point) >= result))
            continue;
         const Int_t left = fBVHChildren[2 * treeNode];
         const Int_t right = fBVHChildren[2 * treeNode + 1];
         if (left >= 0) {
            stack[stackSize++] = right;
            stack[stackSize++] = left;
            continue;
         }
         const Bool_t nodeCanContain = ContainsBox(treeBox, point);
         TransformPointToNode(right, point, local);
         const Bool_t nodeInside = nodeCanContain && NodeContains(right, local);
         if (nodeInside) {
            if (!in)
               return 0.;
            result = std::min(result, NodeSafety(right, local, kTRUE));
         } else if (!in) {
            result = std::min(result, NodeSafety(right, local, kFALSE));
         }
      }
      return result == TGeoShape::Big() ? 0. : result;
   }
   const Bool_t useSmallUnionFastPath = !in && GetNnodes() <= kSmallUnionNodeLimit && fVoxelized &&
                                        fBoxes.size() == static_cast<std::size_t>(kBoxStride * GetNnodes());
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      if (in && !AcceptNode(inode, point))
         continue;
      Bool_t nodeCanContain = kTRUE;
      if (useSmallUnionFastPath) {
         const Double_t *box = &fBoxes[kBoxStride * inode];
         const Double_t tolerance = TGeoShape::Tolerance();
         Double_t boxSafety = 0.;
         for (Int_t axis = 0; axis < 3; ++axis) {
            const Double_t gap = std::max(box[2 * axis] - point[axis], point[axis] - box[2 * axis + 1]);
            if (gap > tolerance) {
               nodeCanContain = kFALSE;
               boxSafety = std::max(boxSafety, gap - tolerance);
            }
         }
         if (boxSafety >= result)
            continue;
      }
      TransformPointToNode(inode, point, local);
      const Bool_t nodeInside = nodeCanContain && NodeContains(inode, local);
      if (nodeInside) {
         if (!in)
            return 0.;
         result = std::min(result, NodeSafety(inode, local, kTRUE));
         continue;
      }
      if (!in)
         result = std::min(result, NodeSafety(inode, local, kFALSE));
   }
   return result == TGeoShape::Big() ? 0. : result;
}

void TGeoMultiUnion::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
{
   for (Int_t i = 0; i < vecsize; ++i)
      safe[i] = Safety(&points[3 * i], inside[i]);
}

Double_t TGeoMultiUnion::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step,
                                         Double_t *safe) const
{
   if (iact < 3 && safe) {
      *safe = Safety(point, kFALSE);
      if (iact == 0)
         return TGeoShape::Big();
      if (iact == 1 && step < *safe)
         return TGeoShape::Big();
   }
   if (TGeoBBox::DistFromOutside(point, dir, fDX, fDY, fDZ, fOrigin, step) >= step)
      return TGeoShape::Big();

   Double_t result = TGeoShape::Big();
   Double_t local[3], localDir[3];
   if (HasBVH()) {
      std::array<Int_t, kBVHStackSize> stack;
      std::size_t stackSize = 1;
      stack[0] = 0;
      while (stackSize) {
         const Int_t treeNode = stack[--stackSize];
         const Double_t limit = std::min(step, result);
         if (BoxRayDistance(&fBVHBoxes[kBoxStride * treeNode], point, dir, limit) >= limit)
            continue;
         const Int_t left = fBVHChildren[2 * treeNode];
         const Int_t right = fBVHChildren[2 * treeNode + 1];
         if (left >= 0) {
            stack[stackSize++] = right;
            stack[stackSize++] = left;
            continue;
         }
         TransformToNode(right, point, dir, local, localDir);
         result = std::min(result, NodeDistFromOutside(right, local, localDir, limit));
      }
      return result < step ? result : TGeoShape::Big();
   }
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      if (!CrossesNodeBox(inode, point, dir, std::min(step, result)))
         continue;
      TransformToNode(inode, point, dir, local, localDir);
      result = std::min(result, NodeDistFromOutside(inode, local, localDir, std::min(step, result)));
   }
   return result < step ? result : TGeoShape::Big();
}

void TGeoMultiUnion::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                                       Double_t *step) const
{
   for (Int_t i = 0; i < vecsize; ++i)
      dists[i] = DistFromOutside(&points[3 * i], &dirs[3 * i], 3, step[i]);
}

Double_t TGeoMultiUnion::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step,
                                        Double_t *safe) const
{
   if (iact < 3 && safe) {
      *safe = Safety(point, kTRUE);
      if (iact == 0)
         return TGeoShape::Big();
      if (iact == 1 && step < *safe)
         return TGeoShape::Big();
   }

   Double_t current[3];
   std::memcpy(current, point, sizeof(current));
   Double_t travelled = 0.;
   constexpr Int_t maxIterations = 10000;
   for (Int_t iteration = 0; iteration < maxIterations; ++iteration) {
      Double_t next = TGeoShape::Big();
      Double_t local[3], localDir[3];
      if (HasBVH()) {
         std::array<Int_t, kBVHStackSize> stack;
         std::size_t stackSize = 1;
         stack[0] = 0;
         while (stackSize) {
            const Int_t treeNode = stack[--stackSize];
            if (!ContainsBox(&fBVHBoxes[kBoxStride * treeNode], current))
               continue;
            const Int_t left = fBVHChildren[2 * treeNode];
            const Int_t right = fBVHChildren[2 * treeNode + 1];
            if (left >= 0) {
               stack[stackSize++] = right;
               stack[stackSize++] = left;
               continue;
            }
            TransformToNode(right, current, dir, local, localDir);
            if (NodeContains(right, local))
               next = std::min(next, NodeDistFromInside(right, local, localDir));
         }
      } else {
         for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
            if (!AcceptNode(inode, current))
               continue;
            TransformToNode(inode, current, dir, local, localDir);
            if (NodeContains(inode, local))
               next = std::min(next, NodeDistFromInside(inode, local, localDir));
         }
      }
      if (next == TGeoShape::Big())
         return travelled;

      const Double_t boundary = travelled + std::max(0., next);
      if (boundary > step)
         return TGeoShape::Big();
      const Double_t push = PushDistance(boundary);
      for (Int_t axis = 0; axis < 3; ++axis)
         current[axis] = point[axis] + (boundary + push) * dir[axis];
      if (!Contains(current))
         return boundary;
      travelled = boundary + push;
   }
   Error("DistFromInside", "Navigation did not converge for multi-union %s", GetName());
   return TGeoShape::Big();
}

void TGeoMultiUnion::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                                      Double_t *step) const
{
   for (Int_t i = 0; i < vecsize; ++i)
      dists[i] = DistFromInside(&points[3 * i], &dirs[3 * i], 3, step[i]);
}

void TGeoMultiUnion::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const
{
   norm[0] = norm[1] = 0.;
   norm[2] = 1.;
   Double_t best = TGeoShape::Big();
   Double_t local[3], localDir[3], localNorm[3], masterNorm[3];
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      TransformToNode(inode, point, dir, local, localDir);
      const Bool_t inside = GetShape(inode)->Contains(local);
      const Double_t distance = GetShape(inode)->Safety(local, inside);
      if (distance > best)
         continue;
      GetShape(inode)->ComputeNormal(local, localDir, localNorm);
      GetMatrix(inode)->LocalToMasterVect(localNorm, masterNorm);
      const Double_t dot = masterNorm[0] * dir[0] + masterNorm[1] * dir[1] + masterNorm[2] * dir[2];
      if (dot < 0.) {
         masterNorm[0] = -masterNorm[0];
         masterNorm[1] = -masterNorm[1];
         masterNorm[2] = -masterNorm[2];
      }
      std::memcpy(norm, masterNorm, sizeof(masterNorm));
      best = distance;
   }
}

void TGeoMultiUnion::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
{
   for (Int_t i = 0; i < vecsize; ++i)
      ComputeNormal(&points[3 * i], &dirs[3 * i], &norms[3 * i]);
}

Double_t TGeoMultiUnion::Capacity() const
{
   if (GetNnodes() == 0)
      return 0.;
   TRandom3 random(0);
   constexpr Int_t samples = 100000;
   Int_t inside = 0;
   Double_t point[3];
   for (Int_t i = 0; i < samples; ++i) {
      point[0] = fOrigin[0] - fDX + 2. * fDX * random.Rndm();
      point[1] = fOrigin[1] - fDY + 2. * fDY * random.Rndm();
      point[2] = fOrigin[2] - fDZ + 2. * fDZ * random.Rndm();
      inside += Contains(point);
   }
   return 8. * fDX * fDY * fDZ * inside / samples;
}

Int_t TGeoMultiUnion::DistancetoPrimitive(Int_t px, Int_t py)
{
   return ShapeDistancetoPrimitive(8, px, py);
}

TGeoVolume *TGeoMultiUnion::Divide(TGeoVolume *, const char *, Int_t, Int_t, Double_t, Double_t)
{
   Error("Divide", "Multi-union shapes cannot be divided");
   return nullptr;
}

void TGeoMultiUnion::InspectShape() const
{
   printf("*** TGeoMultiUnion : %s, %d nodes, voxelized=%s\n", GetName(), GetNnodes(), fVoxelized ? "yes" : "no");
   TGeoBBox::InspectShape();
}

void TGeoMultiUnion::SavePrimitive(std::ostream &out, Option_t *option)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      GetShape(inode)->SavePrimitive(out, option);
      GetMatrix(inode)->SavePrimitive(out, option);
   }
   out << "   TGeoMultiUnion *" << GetPointerName() << " = new TGeoMultiUnion(\"" << GetName() << "\");\n";
   for (Int_t inode = 0; inode < GetNnodes(); ++inode)
      out << "   " << GetPointerName() << "->AddNode(" << GetShape(inode)->GetPointerName() << ", "
          << GetMatrix(inode)->GetPointerName() << ");\n";
   if (fVoxelized)
      out << "   " << GetPointerName() << "->Voxelize();\n";
   TObject::SetBit(kGeoSavePrimitive);
}
