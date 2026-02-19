// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGCone
#define ROOT_TGeoVGCone

// TGeoVGCone
//
// Class description:
//
// Wrapper class for TGeoCone to make use of VecGeom Cone.

#include "TGeoVGAdapter.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class GenericUnplacedCone : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGCone methods
    double GetDz() const { return 0; }
    double GetRmin1() const { return 0; }
    double GetRmax1() const { return 0; }
    double GetRmin2() const { return 0; }
    double GetRmax2() const { return 0; }
    double_t GetPhi1() const { return 0; }
    double_t GetPhi2() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedCone.h>
#endif

class TGeoVGCone final : public TGeoVGAdapter<vecgeom::GenericUnplacedCone> {
                            // in GenericUnplacedCone is used in G$
  using Shape_t = vecgeom::GenericUnplacedCone;
  using Base_t = TGeoVGAdapter<vecgeom::GenericUnplacedCone>;

protected:
   // // data members
   // Double_t fDz;    // half length
   // Double_t fRmin1; // inner radius at -dz
   // Double_t fRmax1; // outer radius at -dz
   // Double_t fRmin2; // inner radius at +dz
   // Double_t fRmax2; // outer radius at +dz
   //                  // methods
   TGeoVGCone(const TGeoVGCone &) = delete;
   TGeoVGCone &operator=(const TGeoVGCone &) = delete;

public:
   // constructors
   TGeoVGCone();
   TGeoVGCone(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2);
   TGeoVGCone(const char *name, Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2);
   TGeoVGCone(Double_t *params);
   // destructor
   ~TGeoVGCone() override;
   // methods

   // static methods
   static Double_t Capacity(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2);
   static void ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t dz, Double_t rmin1,
                              Double_t rmax1, Double_t rmin2, Double_t rmax2);
   static void DistToCone(const Double_t *point, const Double_t *dir, Double_t dz, Double_t r1, Double_t r2,
                          Double_t &b, Double_t &delta);
   static Double_t DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                   Double_t rmax1, Double_t rmin2, Double_t rmax2);
   static Double_t DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                    Double_t rmax1, Double_t rmin2, Double_t rmax2);
   static Double_t SafetyS(const Double_t *point, Bool_t in, Double_t dz, Double_t rmin1, Double_t rmax1,
                           Double_t rmin2, Double_t rmax2, Int_t skipz = 0);

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;

   const char *GetAxisName(Int_t iaxis) const override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   Int_t GetByteCount() const override { return 56; }
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override;

   Double_t GetDz() const { return Base_t::GetDz(); } 
   Double_t GetRmin1() const { return Base_t::GetRmin1(); }
   Double_t GetRmax1() const { return Base_t::GetRmax1(); }
   Double_t GetRmin2() const { return Base_t::GetRmin2(); }
   Double_t GetRmax2() const { return Base_t::GetRmax2(); }

   void InspectShape() const override;
   Bool_t IsCylType() const override { return kTRUE; }
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetConeDimensions(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buffer) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGCone, 1) // conical tube class
};

// TGeoVGConeSeg
//
// Class description:
//
// Wrapper class for TGeoConeSeg to make use of VecGeom Cone.

class TGeoVGConeSeg final : public TGeoVGAdapter<vecgeom::GenericUnplacedCone> {
                            // in GenericUnplacedCone is used in G$
  using Shape_t = vecgeom::GenericUnplacedCone;
  using Base_t = TGeoVGAdapter<vecgeom::GenericUnplacedCone>;

protected:
   // data members
   // Double_t fPhi1; // first phi limit
   // Double_t fPhi2; // second phi limit
   // Transient trigonometric data
   Double_t fS1;   // sin(phi1)
   Double_t fC1;   // cos(phi1)
   Double_t fS2;   // sin(phi2)
   Double_t fC2;   // cos(phi2)
   Double_t fSm;   // sin(0.5*(phi1+phi2))
   Double_t fCm;   // cos(0.5*(phi1+phi2))
   Double_t fCdfi; // cos(0.5*(phi1-phi2))

   void InitTrigonometry();

public:
   // constructors
   TGeoVGConeSeg();
   TGeoVGConeSeg(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t phi1,
               Double_t phi2);
   TGeoVGConeSeg(const char *name, Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2,
               Double_t phi1, Double_t phi2);
   TGeoVGConeSeg(Double_t *params);
   // destructor
   ~TGeoVGConeSeg() override;

   // static methods
   static Double_t
   Capacity(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t phi1, Double_t phi2);
   static void ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t dz, Double_t rmin1,
                              Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1, Double_t s1, Double_t c2,
                              Double_t s2);
   static Double_t DistToCons(const Double_t *point, const Double_t *dir, Double_t r1, Double_t z1, Double_t r2,
                              Double_t z2, Double_t phi1, Double_t phi2);
   static Double_t DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                   Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1, Double_t s1,
                                   Double_t c2, Double_t s2, Double_t cm, Double_t sm, Double_t cdfi);
   static Double_t DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                    Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1, Double_t s1,
                                    Double_t c2, Double_t s2, Double_t cm, Double_t sm, Double_t cdfi);
   static Double_t SafetyS(const Double_t *point, Bool_t in, Double_t dz, Double_t rmin1, Double_t rmax1,
                           Double_t rmin2, Double_t rmax2, Double_t phi1, Double_t phi2, Int_t skipz = 0);

   // methods
   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   // void AfterStreamer() override;
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 64; }
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override;

   Double_t GetDz() const { return Base_t::GetDz(); } 
   Double_t GetRmin1() const { return Base_t::GetRmin1(); }
   Double_t GetRmax1() const { return Base_t::GetRmax1(); }
   Double_t GetRmin2() const { return Base_t::GetRmin2(); }
   Double_t GetRmax2() const { return Base_t::GetRmax2(); }
   Double_t GetPhi1() const;
   Double_t GetPhi2() const;

   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetConsDimensions(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t phi1,
                          Double_t phi2);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buffer) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGConeSeg, 2) // conical tube segment class
};

#endif

#endif
