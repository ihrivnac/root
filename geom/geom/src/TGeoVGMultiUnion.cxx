// @(#)root/geom:$Id$// Author: Andrei Gheata   24/10/01

// Contains() and DistFromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoBBox.h"
#include "TGeoVGMultiUnion.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TGeoBaseBox.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"

#include "TGeoVGArb8.h"
#include "TGeoVGBBox.h"
#include "TGeoVGCone.h"
#include "TGeoVGEltu.h"
#include "TGeoVGHype.h"
#include "TGeoVGPcon.h"
#include "TGeoVGPgon.h"
#include "TGeoVGSphere.h"
#include "TGeoVGTorus.h"
#include "TGeoVGTrd.h"
#include "TGeoVGTube.h"
#include "TGeoVGXtru.h"

// ClassImp(TGeoVGMultiUnion);

namespace {

////////////////////////////////////////////////////////////////////////////////
/// Converts a TGeoHMatrix to a vecgeom::Transformation3D
/// @param matrix The input TGeoMatrix containing rotation and translation.
/// @return A new vecgeom::Transformation3D representing the same transformation.

vecgeom::Transformation3D* ConvertMatrix(const TGeoMatrix& matrix) 
{
   Double_t const *const t = matrix.GetTranslation();
   Double_t const *const r = matrix.GetRotationMatrix();
   vecgeom::cxx::Transformation3D *const transformation =
      new vecgeom::cxx::Transformation3D(t[0], t[1], t[2], r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]);
   return transformation;
}

}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGMultiUnion::TGeoVGMultiUnion()
  : Base_t("")
{
   SetShapeBit(TGeoShape::kGeoComb);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor where half-lengths are provided.

TGeoVGMultiUnion::TGeoVGMultiUnion(const char *name)
  : Base_t(name)
{
   SetShapeBit(TGeoShape::kGeoComb);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVGMultiUnion::~TGeoVGMultiUnion() {}


////////////////////////////////////////////////////////////////////////////////
/// Add component node

void TGeoVGMultiUnion::AddNode(TGeoShape const* shape, const TGeoHMatrix& transform)
{
   if (! shape->IsVecGeom()) {
     Error("AddNode", "shape is not a vecgeom shape");
     std::cerr << "    " << shape->GetName() << std::endl;
     return;
   }

   // Identify the concrete adapted type
   // We use ROOT's Class() mechanism to identify which adapter we have.
   vecgeom::VUnplacedVolume const* unplaced = nullptr;
   if (shape->IsA() == TGeoVGBBox::Class()) {
      // std::cout << "adding VGBBox component" << std::endl;
      auto* adapter = static_cast<TGeoVGBBox const*>(shape);
      // unplaced = static_cast<vecgeom::UnplacedBox const*>(adapter);
      unplaced = adapter->GetUnplaced();
   } 
   else if (shape->IsA() == TGeoVGTube::Class()) {
      auto* adapter = static_cast<TGeoVGTube const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGTubeSeg::Class()) {
      auto* adapter = static_cast<TGeoVGTubeSeg const*>(shape);
       unplaced = adapter->GetUnplaced();
  }
   // to be skipped
   else if (shape->IsA() == TGeoVGCtub::Class()) {
      auto* adapter = static_cast<TGeoVGCtub const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGTorus::Class()) {
      auto* adapter = static_cast<TGeoVGTorus const*>(shape);
      unplaced = adapter->GetUnplaced();
  }
   else if (shape->IsA() == TGeoVGCone::Class()) {
      auto* adapter = static_cast<TGeoVGCone const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   // to be skipped 
   else if (shape->IsA() == TGeoVGConeSeg::Class()) {
      auto* adapter = static_cast<TGeoVGConeSeg const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   // to be skipped
   else if (shape->IsA() == TGeoVGPgon::Class()) {
      auto* adapter = static_cast<TGeoVGPgon const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGPcon::Class()) {
      auto* adapter = static_cast<TGeoVGPcon const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGTrd::Class()) {
      auto* adapter = static_cast<TGeoVGTrd const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGXtru::Class()) {
      auto* adapter = static_cast<TGeoVGXtru const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGArb8::Class()) {
      auto* adapter = static_cast<TGeoVGArb8 const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGTrap::Class()) {
      auto* adapter = static_cast<TGeoVGTrap const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   // to be skipped
   else if (shape->IsA() == TGeoVGGtra::Class()) {
      auto* adapter = static_cast<TGeoVGGtra const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   else if (shape->IsA() == TGeoVGEltu::Class()) {
      auto* adapter = static_cast<TGeoVGEltu const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   // to be skipped
   else if (shape->IsA() == TGeoVGSphere::Class()) {
      auto* adapter = static_cast<TGeoVGSphere const*>(shape);
      unplaced = adapter->GetUnplaced();
   }
   // to be skipped
   else if (shape->IsA() == TGeoVGHype::Class()) {
      auto* adapter = static_cast<TGeoVGHype const*>(shape);
      unplaced = adapter->GetUnplaced();
   }   

   // auto vecGeomShape = dynamic_cast<vecgeom::VUnplacedVolume const*>(shape);
   if (unplaced == nullptr) {
      Error("AddNode", "shape is not supported");
      std::cerr << "    " << shape->GetName() << std::endl;
      return;
   } 

   auto matrix = ConvertMatrix(transform);
 
   Shape_t::AddNode(unplaced, *matrix);

   delete matrix;
}

////////////////////////////////////////////////////////////////////////////////
/// Close VecGeom shape and compute bounding box

void TGeoVGMultiUnion::Close()
{
   Shape_t::Close();
   ComputeBBox();
   SetShapeBit(TGeoShape::kGeoClosedShape);     
}

////////////////////////////////////////////////////////////////////////////////
/// Define Bonding Box via VecGeomShape

void TGeoVGMultiUnion::ComputeBBox()
{
  vecgeom::Vector3D<Double_t> amin, amax;
  // This calls VecGeom to get the extent
  Shape_t::Extent(amin, amax); 

  Double_t dx = (amax[0] - amin[0]) * 0.5;
  Double_t dy = (amax[1] - amin[1]) * 0.5;
  Double_t dz = (amax[2] - amin[2]) * 0.5;
  
  Double_t origin[3];
  origin[0] = (amax[0] + amin[0]) * 0.5;
  origin[1] = (amax[1] + amin[1]) * 0.5;
  origin[2] = (amax[2] + amin[2]) * 0.5;

  // Set both dimensions AND the origin offset
  Base_t::SetBoxDimensions(dx, dy, dz, origin);
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this box shape belonging to volume "voldiv" into ndiv equal volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume. In case a wrong division axis is supplied,
/// returns pointer to volume to be divided.

TGeoVolume *
TGeoVGMultiUnion::Divide(TGeoVolume * /*voldiv*/, const char */*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/, Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Cannot divide VG multiunion");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// FillBuffer3D

void TGeoVGMultiUnion::FillBuffer3D(TBuffer3D &buffer, Int_t reqSections, Bool_t localFrame) const
{
   // 1. Fill Core and Bounding Box sections using the base TGeoBBox logic
   // This handles the name, color, and the fDX, fDY, fDZ you set in ComputeBBox
   TGeoShape::FillBuffer3D(buffer, reqSections, localFrame);

   // 2. If the viewer wants the actual mesh (the shape's faces)
   if (reqSections & TBuffer3D::kRaw) {
      // For a MultiUnion, providing a full raw mesh manually is complex.
      // Strategy: If you have many nodes, ROOT prefers that you 
      // handle this via a TGeoCompositeShape or by filling 
      // the buffer with the vertices of the bounding box as a placeholder.
      
      buffer.SetSectionsValid(TBuffer3D::kRaw);
      // Note: To see the individual parts, you would ideally iterate 
      // through the VecGeom nodes and merge their polyhedra.
   }
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGMultiUnion::GetMakeRuntimeShape(TGeoShape */*mother*/, TGeoMatrix */*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;

   Error("GetMakeRuntimeShape", "Cannot make runtime shape for VG multiunion");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGMultiUnion::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const 
{
// #ifndef __ROOTCLING__
   size_t numNodes = Shape_t::GetNumberOfSolids();
   // Estimation: sum of vertices of components
   // A safe upper bound for boxes is 8 per node.
   nvert = numNodes * 8; 
   nsegs = numNodes * 12;
   npols = numNodes * 6;
// #else
//    nvert = 8; nsegs = 12; npols = 6;
// #endif
}

////////////////////////////////////////////////////////////////////////////////
/// Returns NmeshVertices

Int_t TGeoVGMultiUnion::GetNmeshVertices() const 
{
// #ifndef __ROOTCLING__
   // Just return a safe estimate based on components to avoid complex mesh logic
   return Shape_t::GetNumberOfSolids() * 8;
// #else
//    return 8;
// #endif   
}

////////////////////////////////////////////////////////////////////////////////
/// Prints shape parameters

void TGeoVGMultiUnion::InspectShape() const
{
   printf("*** Shape %s: TGeoVGMultiUnion ***\n", GetName());
   
   // 1. Print Bounding Box Info (calculated in ComputeBBox)
   printf("    Bounding Box: DX=%7.3f DY=%7.3f DZ=%7.3f\n", 
          fBoundingBox.GetDX(), fBoundingBox.GetDY(), fBoundingBox.GetDZ());
   printf("    Origin:       X=%7.3f Y=%7.3f Z=%7.3f\n", 
          fBoundingBox.GetOrigin()[0], fBoundingBox.GetOrigin()[1], fBoundingBox.GetOrigin()[2]);

//#ifndef __ROOTCLING__
   // 2. Access the VecGeom internal nodes
   // Get the number of nodes added to the MultiUnion
   size_t numNodes = Shape_t::GetNumberOfSolids(); 
   printf("    Number of components: %zu\n", numNodes);

   for (size_t i = 0; i < numNodes; ++i) {
      // Retrieve the placed volume (node)
      vecgeom::VPlacedVolume const* placed = Shape_t::GetNode(i);
      vecgeom::VUnplacedVolume const* unplaced = placed->GetUnplacedVolume();
      
      // Get Transformation details
      vecgeom::Transformation3D const* transform = placed->GetTransformation();
      vecgeom::Vector3D<vecgeom::Precision> transl = transform->Translation();

      printf("    - Node [%2zu]: \n", i);
      printf("             ShapeType: %d \n", static_cast<int>(unplaced->GetType()));
      printf("             Translation: (%7.3f, %7.3f, %7.3f)\n", 
             transl.x(), transl.y(), transl.z());
      
      // Optional: Ask the specific unplaced volume to print itself
      // unplaced->Print(); 
   }
// #else
//    printf("    [Details hidden for ROOT dictionary generation]\n");
// #endif
   printf("**************************************************\n");
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGMultiUnion::MakeBuffer3D() const
{
   // Create a standard buffer for a composite-like shape
   TBuffer3D *buff = new TBuffer3D(TBuffer3DTypes::kComposite, 8, 24, 12, 0, 0, 0);
   if (buff) FillBuffer3D(*buff, TBuffer3D::kCore | TBuffer3D::kBoundingBox, kTRUE);
   return buff;
}

///
////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGMultiUnion::SavePrimitive(std::ostream &/*out*/, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;

   Error("SavePrimitive", "Not implemented fir VG Multiunion");
}

////////////////////////////////////////////////////////////////////////////////
/// Fill box vertices to an array.

void TGeoVGMultiUnion::SetBoxPoints(Double_t *points) const
{
   Base_t::SetPoints(points);
}

#endif
