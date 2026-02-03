// @(#)root/base:$Id$

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGTubeSeg
#define ROOT_TGeoVGTubeSeg

// TGeoVGTubeSeg
//
// Class description:
//
// Wrapper class for TGeoVGTubeSeg to make use of VecGeom Tube.

#include "TGeoVGAdapter.h"

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class GenericUnplacedTube : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGTube methods
    double rmin() const { return 0; }
    double rmax() const { return 0; }
    double z() const { return 0; }
    double sphi() const { return 0; }
    double dphi() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedTube.h>
#endif

class TGeoVGTubeSeg : public TGeoVGAdapter<vecgeom::GenericUnplacedTube> {  
                            // in GenericUnplacedTube is used in G$
  using Shape_t = vecgeom::GenericUnplacedTube;
  using Base_t = TGeoVGAdapter<vecgeom::GenericUnplacedTube>;

protected:
   // data members
   // Transient trigonometric data (use to compute BBox)
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
   TGeoVGTubeSeg();
   TGeoVGTubeSeg(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2);
   TGeoVGTubeSeg(const char *name, Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2);
   TGeoVGTubeSeg(Double_t *params);

   // destructor
   ~TGeoVGTubeSeg() {}

   // static methods
   static Double_t Capacity(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2);
   static void ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t rmin, Double_t rmax,
                              Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2);
   static Double_t DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax,
                                   Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2, Double_t cm,
                                   Double_t sm, Double_t cdfi);
   static Double_t DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax,
                                    Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2, Double_t cm,
                                    Double_t sm, Double_t cdfi);
   static Double_t SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1,
                           Double_t phi2, Int_t skipz = 0);

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 56; }
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   Int_t GetNmeshVertices() const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;

   Double_t GetRmin() const { return rmin(); }
   Double_t GetRmax() const { return rmax(); }
   Double_t GetDz() const { return z(); }
   Double_t GetPhi1() const;
   Double_t GetPhi2() const;
   void SetTubsDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGTubeSeg, 2) // cylindrical tube segment class
};

#endif
