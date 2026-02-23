// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGSphere
#define ROOT_TGeoVGSphere

// TGeoVGSphere
//
// Class description:
//
// Wrapper class for TGeoSphere to make use of VecGeom Shere.

#include "TGeoVGAdapter.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedSphere : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGSphere methods
    double GetInnerRadius() const { return 0; }
    double GetOuterRadius() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedSphere.h>
#endif

class TGeoVGSphere final : public TGeoVGAdapter<vecgeom::UnplacedSphere> {
                            // in UnplacedSphere is used in G$
  using Shape_t = vecgeom::UnplacedSphere;
  using Base_t = TGeoVGAdapter<vecgeom::UnplacedSphere>;

protected:
   Int_t fNz;        // number of z planes for drawing
   Int_t fNseg;      // number of segments for drawing

   TGeoVGSphere(const TGeoVGSphere &) = delete;
   TGeoVGSphere &operator=(const TGeoVGSphere &) = delete;

public:
   // constructors
   TGeoVGSphere();
   TGeoVGSphere(Double_t rmin, Double_t rmax, Double_t theta1 = 0, Double_t theta2 = 180, Double_t phi1 = 0,
              Double_t phi2 = 360);
   TGeoVGSphere(const char *name, Double_t rmin, Double_t rmax, Double_t theta1 = 0, Double_t theta2 = 180,
              Double_t phi1 = 0, Double_t phi2 = 360);
   TGeoVGSphere(Double_t *param, Int_t nparam = 6);
   // destructor
   ~TGeoVGSphere() override;

   // methods
   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   // Double_t DistToSphere(const Double_t *point, const Double_t *dir, Double_t rsph, Bool_t check = kTRUE,
   //                       Bool_t firstcross = kTRUE) const;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   const char *GetAxisName(Int_t iaxis) const override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 42; }
   TGeoShape *GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const override { return nullptr; }
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   Int_t GetNumberOfDivisions() const { return fNseg; }
   Bool_t GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const override { return kFALSE; }

   Int_t GetNz() const { return fNz; }
   virtual Double_t GetRmin() const { return GetInnerRadius(); }
   virtual Double_t GetRmax() const { return GetOuterRadius(); }
   Double_t GetTheta1() const;
   Double_t GetTheta2() const;
   Double_t GetPhi1() const;
   Double_t GetPhi2() const;

   void InspectShape() const override;
   Bool_t IsCylType() const override { return kFALSE; }
   Int_t IsOnBoundary(const Double_t *point) const;
   Bool_t
   IsPointInside(const Double_t *point, Bool_t checkR = kTRUE, Bool_t checkTh = kTRUE, Bool_t checkPh = kTRUE) const;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetSphDimensions(Double_t rmin, Double_t rmax, Double_t theta1, Double_t theta2, Double_t phi1, Double_t phi2);
   virtual void SetNumberOfDivisions(Int_t p);
   void SetDimensions(Double_t *param) override;
   void SetDimensions(Double_t *param, Int_t nparam);
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGSphere, 1) // sphere class
};

#endif

#endif
