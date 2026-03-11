// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGMultiUnion
#define ROOT_TGeoVGMultiUnion

// 
//
// Class description:
//
// Wrapper class for TGeoVGMultiUnion to make use of VecGeom MultiUnion.

#include "TGeoVGAdapter.h"

namespace vecgeom {
inline namespace cxx {
class Transformation3D;
}
}

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedMultiUnion : public VUnplacedVolume {
  };
}
#else
#include <VecGeom/volumes/UnplacedMultiUnion.h>
#endif

class TGeoVGMultiUnion final : public TGeoVGAdapter<vecgeom::UnplacedMultiUnion>
{
  using Shape_t = vecgeom::UnplacedMultiUnion;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedMultiUnion>;

protected:
   void FillBuffer3D(TBuffer3D &buffer, Int_t reqSections, Bool_t localFrame) const override;

   TGeoVGMultiUnion(const TGeoVGMultiUnion &) = delete;
   TGeoVGMultiUnion &operator=(const TGeoVGMultiUnion &) = delete;

public:
   // constructors
   TGeoVGMultiUnion();
   TGeoVGMultiUnion(const char *name);
   // destructor
   ~TGeoVGMultiUnion() override;
   
   // Method for solid construction
#ifndef __ROOTCLING__
   using Shape_t::AddNode;
   using Shape_t::Close;
#endif

   void AddNode(TGeoShape const* shape, const TGeoHMatrix& transform);
   void Close();

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;

   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetBoxPoints(Double_t *points) const;

   // ClassDefOverride(TGeoVGMultiUnion, 1) // box primitive
};

#endif

#endif
