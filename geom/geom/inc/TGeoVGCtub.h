// @(#)root/base:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGCtub
#define ROOT_TGeoVGCtub

// TGeoVGCtub
//
// Class description:
//
// Wrapper class for TGeoCtub to make use of VecGeom Tube.

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
#include <VecGeom/volumes/UnplacedCutTube.h>
#endif

class TGeoVGCtub : public TGeoVGAdapter<vecgeom::UnplacedCutTube> {
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
   // Arrays for normals
   Double_t fNlow[3];  // normal to lower cut plane
   Double_t fNhigh[3]; // normal to higher cut plane

   void InitTrigonometry();  // from TgeoTubeSeg

public:
   // constructors
   TGeoVGCtub();
   TGeoVGCtub(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2, Double_t lx, Double_t ly,
            Double_t lz, Double_t tx, Double_t ty, Double_t tz);
   TGeoVGCtub(const char *name, Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2, Double_t lx,
            Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz);
   TGeoVGCtub(Double_t *params);
   // destructor
   ~TGeoVGCtub() override;

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   // methods
   // Double_t Capacity() const override;
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 98; }
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   void InspectShape() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;

   Double_t GetRmin() const { return rmin(); }
   Double_t GetRmax() const { return rmax(); }
   Double_t GetDz() const { return z(); }
   Double_t GetPhi1() const { return sphi(); }
   Double_t GetPhi2() const { return sphi() + dphi(); }
   const Double_t *GetNlow() const { return &fNlow[0]; }
   const Double_t *GetNhigh() const { return &fNhigh[0]; }
   Double_t GetZcoord(Double_t xc, Double_t yc, Double_t zc) const;

   void SetTubsDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2);
   void SetCtubDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2, Double_t lx,
                          Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz);
   void SetDimensions(Double_t *param);
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;

   ClassDefOverride(TGeoVGCtub, 1) // cut tube segment class
};

#endif
