// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGMultiUnion
#define ROOT_TGeoVGMultiUnion

// TGeoVGPcon
//
// Class description:
//
// Wrapper class for TGeoPcon to make use of VecGeom Pcon.

#include "TGeoVGAdapter.h"

#include <iostream>
#include <vector>

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedMultiUnion : public VUnplacedVolume {
    // nothing to be defined    
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedMultiUnion.h>
#endif

class TGeoShape;
class TGeoHMatrix;
class TGeoCompositeShape;

namespace vecgeom {
inline namespace cxx {
class Transformation3D;
}
}

typedef std::pair<TGeoShape*, TGeoHMatrix*> ShapeEntry;

class TGeoVGMultiUnion final : public TGeoVGAdapter<vecgeom::UnplacedMultiUnion> {
  using Shape_t = vecgeom::UnplacedMultiUnion;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedMultiUnion>;

private:
   // data members
   TGeoCompositeShape *fCompositeShape; // the shape define in old style for other than navigation function

protected:
   TGeoVGMultiUnion(const TGeoVGMultiUnion &) = delete;
   TGeoVGMultiUnion &operator=(const TGeoVGMultiUnion &) = delete;

public:
   // constructors
   TGeoVGMultiUnion();
   TGeoVGMultiUnion(const char *name, const std::vector<ShapeEntry>& components);
   // destructor
   ~TGeoVGMultiUnion() override;

   // static functions
   // static vecgeom::cxx::Transformation3D *Convert(TGeoMatrix const *const geomatrix);

   // methods
   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   // void ClearThreadData() const override;
   // void CreateThreadData(Int_t nthreads) override;
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;

   void GetBoundingCylinder(Double_t * /*param*/) const override {}
   TGeoShape *GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const override { return nullptr; }
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const override { return kFALSE; }
   void InspectShape() const override;
   Bool_t IsComposite() const override { return kTRUE; }
   Bool_t IsCylType() const override { return kFALSE; }
   // virtual Bool_t PaintComposite(Option_t *option = "") const;
   void RegisterYourself();
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetDimensions(Double_t * /*param*/) override {}
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGMultiUnion, 1) // boolean composite shape
};

#endif

#endif
