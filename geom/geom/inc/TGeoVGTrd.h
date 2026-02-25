// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGTrd
#define ROOT_TGeoVGTrd

// TGeoVGTrd
//
// Class description:
//
// Wrapper class for TGeoTrd to make use of VecGeom Trd.

#include "TGeoVGAdapter.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class GenericUnplacedTrd : public VUnplacedVolume {
  public:
    // Define the specific getters used in the TGeoVGTube methods
    double dx1() const { return 0; }
    double dx2() const { return 0; }
    double dy1() const { return 0; }
    double dy2() const { return 0; }
    double dz() const { return 0; }
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedTrd.h>
#endif

class TGeoVGTrd : public TGeoVGAdapter<vecgeom::GenericUnplacedTrd> 
{
  using Shape_t = vecgeom::GenericUnplacedTrd;
  using Base_t  = TGeoVGAdapter<vecgeom::GenericUnplacedTrd>;

   // methods
   TGeoVGTrd(const TGeoVGTrd &) = delete;
   TGeoVGTrd &operator=(const TGeoVGTrd &) = delete;

public:
   // constructors
   TGeoVGTrd();
   // Trd1
   TGeoVGTrd(Double_t dx1, Double_t dx2, Double_t dy, Double_t dz);
   TGeoVGTrd(const char *name, Double_t dx1, Double_t dx2, Double_t dy, Double_t dz);
   // Trd2
   TGeoVGTrd(Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2, Double_t dz);
   TGeoVGTrd(const char *name, Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2, Double_t dz);
   // TGeoVGTrd(Double_t *params); 
      // disabled ctor as it to distinguish Trd1 with 4 parameters and Trd2 with 5 parameters
   // destructor
   ~TGeoVGTrd() override;

   // methods

   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override;
   void GetBoundingCylinder(Double_t *param) const override;
   Int_t GetByteCount() const override { return 56; }

   Double_t GetDx1() const { return Base_t::dx1(); }
   Double_t GetDx2() const { return Base_t::dx2(); }
   Double_t GetDy1() const { return Base_t::dy1(); }
   Double_t GetDy2() const { return Base_t::dy2(); }
   Double_t GetDz() const { return  Base_t::dz(); }
   Double_t GetDy() const;  // returns pnly if Trd1

   Int_t
   GetFittingBox(const TGeoBaseBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   // void GetVisibleCorner(const Double_t *point, Double_t *vertex, Double_t *normals) const;
   // void GetOppositeCorner(const Double_t *point, Int_t inorm, Double_t *vertex, Double_t *normals) const;
   void InspectShape() const override;
   Bool_t IsCylType() const override { return kFALSE; }
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetVertex(Double_t *vertex) const;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGTrd, 1) // TRD2 shape class
};

#endif // ROOT_USE_VECGEOM_SOLIDS

#endif
