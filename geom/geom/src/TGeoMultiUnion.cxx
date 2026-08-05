// Author: ROOT team

#include "TGeoMultiUnion.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TRandom3.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

ClassImp(TGeoMultiUnion);

namespace {
constexpr Int_t kBoxStride = 6;
constexpr Int_t kSmallUnionNodeLimit = 9;

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
   fBoxes.clear();
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
   const auto *matrix = GetMatrix(inode);
   matrix->MasterToLocal(point, local);
   if (dir)
      matrix->MasterToLocalVect(dir, localDir);
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

Bool_t TGeoMultiUnion::CrossesNodeBox(Int_t inode, const Double_t *point, const Double_t *dir, Double_t step) const
{
   if (!fVoxelized || fBoxes.size() != static_cast<std::size_t>(kBoxStride * GetNnodes()))
      return kTRUE;

   const Double_t *box = &fBoxes[kBoxStride * inode];
   Double_t enter = 0.;
   Double_t leave = step;
   for (Int_t axis = 0; axis < 3; ++axis) {
      const Double_t lo = box[2 * axis];
      const Double_t hi = box[2 * axis + 1];
      if (std::abs(dir[axis]) < std::numeric_limits<Double_t>::epsilon()) {
         if (point[axis] < lo || point[axis] > hi)
            return kFALSE;
         continue;
      }
      Double_t first = (lo - point[axis]) / dir[axis];
      Double_t last = (hi - point[axis]) / dir[axis];
      if (first > last)
         std::swap(first, last);
      enter = std::max(enter, first);
      leave = std::min(leave, last);
      if (enter > leave)
         return kFALSE;
   }
   return leave >= 0.;
}

void TGeoMultiUnion::Voxelize()
{
   fBoxes.assign(kBoxStride * GetNnodes(), 0.);
   Double_t vertices[24];
   Double_t master[3];
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      auto *shape = GetShape(inode);
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
   ComputeBBox();
}

void TGeoMultiUnion::ComputeBBox()
{
   if (GetNnodes() == 0) {
      SetBoxDimensions(0., 0., 0.);
      return;
   }

   const Bool_t restoreVoxelized = fVoxelized;
   if (fBoxes.size() != static_cast<std::size_t>(kBoxStride * GetNnodes())) {
      fVoxelized = kFALSE;
      Voxelize();
      fVoxelized = restoreVoxelized;
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
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      if (!AcceptNode(inode, point))
         continue;
      GetMatrix(inode)->MasterToLocal(point, local);
      if (GetShape(inode)->Contains(local))
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
      GetMatrix(inode)->MasterToLocal(point, local);
      const Bool_t nodeInside = nodeCanContain && GetShape(inode)->Contains(local);
      if (nodeInside) {
         if (!in)
            return 0.;
         result = std::min(result, GetShape(inode)->Safety(local, kTRUE));
         continue;
      }
      if (!in)
         result = std::min(result, GetShape(inode)->Safety(local, kFALSE));
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
   for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
      if (!CrossesNodeBox(inode, point, dir, std::min(step, result)))
         continue;
      TransformToNode(inode, point, dir, local, localDir);
      result = std::min(result, GetShape(inode)->DistFromOutside(local, localDir, 3, std::min(step, result)));
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
      for (Int_t inode = 0; inode < GetNnodes(); ++inode) {
         if (!AcceptNode(inode, current))
            continue;
         TransformToNode(inode, current, dir, local, localDir);
         if (GetShape(inode)->Contains(local))
            next = std::min(next, GetShape(inode)->DistFromInside(local, localDir, 3));
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
