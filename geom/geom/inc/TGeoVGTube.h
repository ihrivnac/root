// @(#)root/base:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGTube
#define ROOT_TGeoVGTube

// TGeoVGTube
//
// Class description:
//
// Wrapper class for TGeoTube to make use of VecGeom Tube.

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

class TGeoVGTube : public TGeoVGAdapter<vecgeom::GenericUnplacedTube> {  
                            // in GenericUnplacedTube is used in G$
  using Shape_t = vecgeom::GenericUnplacedTube;
  using Base_t = TGeoVGAdapter<vecgeom::GenericUnplacedTube>;
public:
   // constructors
   TGeoVGTube();
   TGeoVGTube(Double_t rmin, Double_t rmax, Double_t dz);
   TGeoVGTube(const char *name, Double_t rmin, Double_t rmax, Double_t dz);
   TGeoVGTube(Double_t *params);
   TGeoVGTube(const TGeoVGTube &) = delete;
   TGeoVGTube &operator=(const TGeoVGTube &) = delete;

   // destructor
   virtual ~TGeoVGTube() {}
   
   // static methods
   static Double_t Capacity(Double_t rmin, Double_t rmax, Double_t dz);   
   static void ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t rmin, Double_t rmax,
                              Double_t dz);
   static Double_t
   DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz);
   static Double_t
   DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz);
   static void DistToTube(Double_t rsq, Double_t nsq, Double_t rdotn, Double_t radius, Double_t &b, Double_t &delta);
   static Double_t
   SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz, Int_t skipz = 0);

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   const char *GetAxisName(Int_t iaxis) const override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;

   inline Double_t GetRmin() const { return rmin(); }
   inline Double_t GetRmax() const { return rmax(); }
   inline Double_t GetDz() const { return z(); }
   inline Bool_t HasRmin() const { return (rmin() > 0) ? kTRUE : kFALSE; }
   void SetTubeDimensions(Double_t rmin, Double_t rmax, Double_t dz);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGTube, 1) // cylindrical l tube class
};

#endif
