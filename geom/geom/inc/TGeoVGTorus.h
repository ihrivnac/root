// @(#)root/base:$Id$
// Author: Andrei Gheata   28/07/03

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGTorus
#define ROOT_TGeoVGTorus

// TGeoVGTorus
//
// Class description:
//
// Wrapper class for TGeoTorus to make use of VecGeom Torus.

#include "TGeoVGAdapter.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedTorus2 : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGTorus methods
    double rmin() const { return 0; }
    double rmax() const { return 0; }
    double rtor() const { return 0; }
    double sphi() const { return 0; }
    double dphi() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedTorus2.h>
#endif

class TGeoVGTorus final : public TGeoVGAdapter<vecgeom::UnplacedTorus2> {
  using Shape_t = vecgeom::UnplacedTorus2;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedTorus2>;

public:
   // Double_t Capacity() const override;
   Double_t Daxis(const Double_t *pt, const Double_t *dir, Double_t t) const;
   Double_t DDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const;
   Double_t DDDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const;
   Double_t ToBoundary(const Double_t *pt, const Double_t *dir, Double_t r, Bool_t in) const;
   Int_t SolveCubic(Double_t a, Double_t b, Double_t c, Double_t *x) const;
   Int_t SolveQuartic(Double_t a, Double_t b, Double_t c, Double_t d, Double_t *x) const;

public:
   // constructors
   TGeoVGTorus();
   TGeoVGTorus(Double_t r, Double_t rmin, Double_t rmax, Double_t phi1 = 0, Double_t dphi = 360);
   TGeoVGTorus(const char *name, Double_t r, Double_t rmin, Double_t rmax, Double_t phi1 = 0, Double_t dphi = 360);
   TGeoVGTorus(Double_t *params);
   TGeoVGTorus(const TGeoVGTorus &) = delete;
   TGeoVGTorus &operator=(const TGeoVGTorus &) = delete;

   // destructor
   ~TGeoVGTorus() override {}
   // methods

   void ComputeBBox() override;

   // void ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const override;
   // void ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize) override;
   // Bool_t Contains(const Double_t *point) const override;
   // void Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const override;
   // Double_t DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact = 1, Double_t step = TGeoShape::Big(),
   //                         Double_t *safe = nullptr) const override;
   // void DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
   //                       Double_t *step) const override;
   // Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
   //                          Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   // void DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
   //                        Double_t *step) const override;
   // Int_t DistancetoPrimitive(Int_t px, Int_t py) override;

   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   const char *GetAxisName(Int_t iaxis) const override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 56; }
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const override { return kFALSE; }
   void InspectShape() const override;
   Bool_t IsCylType() const override { return kTRUE; }
   TBuffer3D *MakeBuffer3D() const override;

   Double_t GetR() const { return Base_t::rtor(); }
   Double_t GetRmin() const { return Base_t::rmin(); }
   Double_t GetRmax() const { return Base_t::rmax(); }
   Double_t GetPhi1() const;
   Double_t GetDphi() const;

   // Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override;
   // void Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetTorusDimensions(Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGTorus, 1) // torus class
};

#endif // ROOT_USE_VECGEOM_SOLIDS

#endif
