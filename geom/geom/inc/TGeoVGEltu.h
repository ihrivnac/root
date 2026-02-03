// @(#)root/geom:$Id$
// Author: Mihaela Gheata   05/06/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGEltu
#define ROOT_TGeoVGEltu

// TGeoVGEltu
//
// Class description:
//
// Wrapper class for TGeoEltu to make use of VecGeom Eltu.

#include "TGeoVGAdapter.h"

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedEllipticalTube : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGTube methods
    double GetDx() const { return 0; }
    double GetDy() const { return 0; }
    double GetDz() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedEllipticalTube.h>
#endif

class TGeoVGEltu : public TGeoVGAdapter<vecgeom::UnplacedEllipticalTube>
{
  using Shape_t = vecgeom::UnplacedEllipticalTube;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedEllipticalTube>;

public:
   // constructors
   // TGeoVGEltu();
   TGeoVGEltu(Double_t a, Double_t b, Double_t dz);
   TGeoVGEltu(const char *name, Double_t a, Double_t b, Double_t dz);
   TGeoVGEltu(Double_t *params);
   TGeoVGEltu(const TGeoVGEltu &) = delete;
   TGeoVGEltu &operator=(const TGeoVGEltu &) = delete;
   // destructor
   ~TGeoVGEltu() override;

   // methods

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   // Double_t Capacity() const override;
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;

   Double_t GetA() const { return Base_t::GetDx(); }  // TODO: CHECK
   Double_t GetB() const { return Base_t::GetDy(); }
   inline Double_t GetDz() const { return Base_t::GetDz(); }
   Bool_t HasRmin() const { return (Base_t::GetDx() > 0) ? kTRUE : kFALSE; }

   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const override { return kFALSE; }
   void InspectShape() const override;
   Bool_t IsCylType() const override { return kTRUE; }
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetEltuDimensions(Double_t a, Double_t b, Double_t dz);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;

   // ClassDefOverride(TGeoVGEltu, 1) // elliptical tube class
};

#endif
