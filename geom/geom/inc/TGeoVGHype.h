// @(#)root/geom:$Id$
// Author: Mihaela Gheata   20/11/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGHype
#define ROOT_TGeoVGHype

// TGeoVGHype
//
// Class description:
//
// Wrapper class for TGeoHype to make use of VecGeom Hype.

#include "TGeoVGAdapter.h"

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class GenericUnplacedHype : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGTube methods
    double GetRmin() const { return 0.; }
    double GetRmax() const { return 0.; }
    double GetDz() const { return 0.; }
    double GetStIn() const { return 0.; }
    double GetStOut() const { return 0.; }
  };
}
#else
#include <VecGeom/volumes/UnplacedHype.h>
#endif

class TGeoVGHype : public TGeoVGAdapter<vecgeom::GenericUnplacedHype>
{
  using Shape_t = vecgeom::GenericUnplacedHype;
  using Base_t  = TGeoVGAdapter<vecgeom::GenericUnplacedHype>;

private:
   // Precomputed parameters:
   Double_t fTin;    // Tangent of stereo angle for inner surface
   Double_t fTout;   // Tangent of stereo angle for outer surface
   Double_t fTinsq;  // Squared tangent of stereo angle for inner surface
   Double_t fToutsq; // Squared tangent of stereo angle for outer surface

public:
   // constructors
   TGeoVGHype();
   TGeoVGHype(Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz);
   TGeoVGHype(const char *name, Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz);
   TGeoVGHype(Double_t *params);
   TGeoVGHype(const TGeoVGHype &) = delete;
   TGeoVGHype &operator=(const TGeoVGHype &) = delete;
   // destructor
   ~TGeoVGHype() override;
   // methods

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   // Double_t Capacity() const override;
   void ComputeBBox() override;
   // Int_t DistToHype(const Double_t *point, const Double_t *dir, Double_t *s, Bool_t inner, Bool_t in) const;

   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   // Int_t GetByteCount() const override { return 64; }
   Bool_t GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const override { return kFALSE; }
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;

   Double_t GetRmin() const { return Base_t::GetRmin(); }
   Double_t GetRmax() const { return Base_t::GetRmax(); }
   Double_t GetDz() const { return Base_t::GetDz(); }
   Double_t GetStIn() const;
   Double_t GetStOut() const;
   Bool_t HasInner() const { return !TestShapeBit(kGeoRSeg); }
   Bool_t HasRmin() const { return (Base_t::GetRmin() > 0) ? kTRUE : kFALSE; }
   Double_t RadiusHypeSq(Double_t z, Bool_t inner) const;
   Double_t ZHypeSq(Double_t r, Bool_t inner) const;

   void InspectShape() const override;
   Bool_t IsCylType() const override { return kTRUE; }
   TBuffer3D *MakeBuffer3D() const override;
   // virtual void          Paint(Option_t *option);
   // Double_t SafetyToHype(const Double_t *point, Bool_t inner, Bool_t in) const;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetHypeDimensions();
   // void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGHype, 1) // hyperboloid class
};

#endif
