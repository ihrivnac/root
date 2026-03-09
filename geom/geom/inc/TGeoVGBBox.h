// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGBBox
#define ROOT_TGeoVGBBox

// 
//
// Class description:
//
// Wrapper class for TGeoVGBBox to make use of VecGeom Box.

#include "TGeoVGAdapter.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedBox : public VUnplacedVolume {
    double x() const { return 0.; }
    double y() const { return 0.; }
    double z() const { return 0.; }
  };
}
#else
#include <VecGeom/volumes/UnplacedBox.h>
#endif

class TGeoVGBBox final : public TGeoVGAdapter<vecgeom::UnplacedBox>
{
  using Shape_t = vecgeom::UnplacedBox;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedBox>;

protected:
   // void FillBuffer3D(TBuffer3D &buffer, Int_t reqSections, Bool_t localFrame) const override;

   TGeoVGBBox(const TGeoVGBBox &) = delete;
   TGeoVGBBox &operator=(const TGeoVGBBox &) = delete;

public:
   // constructors
   TGeoVGBBox();
   TGeoVGBBox(Double_t dx, Double_t dy, Double_t dz);
   TGeoVGBBox(const char *name, Double_t dx, Double_t dy, Double_t dz);
     // VecGeom does not support shifted origin
   TGeoVGBBox(Double_t *param);
   // destructor
   ~TGeoVGBBox() override;
   // methods
   
   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   static Bool_t Contains(const Double_t *point, Double_t dx, Double_t dy, Double_t dz, const Double_t *origin);
   static Double_t DistFromInside(const Double_t *point, const Double_t *dir, Double_t dx, Double_t dy, Double_t dz,
                                  const Double_t *origin, Double_t stepmax = TGeoShape::Big());
   static Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Double_t dx, Double_t dy, Double_t dz,
                                   const Double_t *origin, Double_t stepmax = TGeoShape::Big());
   // Bool_t CouldBeCrossed(const Double_t *point, const Double_t *dir) const override;

   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override;
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override { return 8; }

   // virtual Double_t GetDX() const { return x(); }
   // virtual Double_t GetDY() const { return y(); }
   // virtual Double_t GetDZ() const { return z(); }
   // virtual const Double_t *GetOrigin() const { return fOrigin; }
   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetBoxDimensions(Double_t dx, Double_t dy, Double_t dz);
   void SetDimensions(Double_t *param) override;
   void SetBoxPoints(Double_t *points) const;

   // ClassDefOverride(TGeoVGBBox, 1) // box primitive
};

#endif

#endif
