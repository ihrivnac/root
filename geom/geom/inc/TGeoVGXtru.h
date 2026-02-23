// @(#)root/geom:$Id$
// Author: Mihaela Gheata   24/01/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGXtru
#define ROOT_TGeoVGXtru

// TGeoVGXtru
//
// Class description:
//
// Wrapper class for TGeoXtru to make use of VecGeom Extruded solid.

#include "TGeoVGAdapter.h"

#include <mutex>
#include <vector>

class TGeoPolygon;

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedExtruded : public VUnplacedVolume {
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedExtruded.h>
#endif

class TGeoVGXtru final : public TGeoVGAdapter<vecgeom::UnplacedExtruded> {
                            // in UnplacedExtruded is used in G$
  using Shape_t = vecgeom::UnplacedExtruded;
  using Base_t = TGeoVGAdapter<vecgeom::UnplacedExtruded>;

public:
   struct ThreadData_t {
      Int_t fSeg;         // !current segment [0,fNvert-1]
      Int_t fIz;          // !current z plane [0,fNz-1]
      Double_t *fXc;      // ![fNvert] current X positions for polygon vertices
      Double_t *fYc;      // ![fNvert] current Y positions for polygon vertices
      TGeoPolygon *fPoly; // !polygon defining section shape

      ThreadData_t();
      ~ThreadData_t();
   };
   ThreadData_t &GetThreadData() const;
   void ClearThreadData() const override;
   void CreateThreadData(Int_t nthreads) override;

protected:
   // data members
    // Data for construction
   Int_t fNvert;       // number of vertices of the 2D polygon (at least 3)
   Int_t fNz;          // number of z planes (at least two)
   Double_t fZcurrent; // current Z position
   Double_t *fX;       //[fNvert] X positions for polygon vertices
   Double_t *fY;       //[fNvert] Y positions for polygon vertices
   Double_t *fZ;       //[fNz] array of Z planes positions
   Double_t *fScale;   //[fNz] array of scale factors (for each Z)
   Double_t *fX0;      //[fNz] array of X offsets (for each Z)
   Double_t *fY0;      //[fNz] array of Y offsets (for each Z)

   mutable std::vector<ThreadData_t *> fThreadData; //! Navigation data per thread
   mutable Int_t fThreadSize;                       //! size of thread-specific array
   mutable std::mutex fMutex;                       //! mutex for thread data

   TGeoVGXtru(const TGeoVGXtru &) = delete;
   TGeoVGXtru &operator=(const TGeoVGXtru &) = delete;

   // methods
   // void GetPlaneVertices(Int_t iz, Int_t ivert, Double_t *vert) const;
   // void GetPlaneNormal(const Double_t *vert, Double_t *norm) const;
   // Bool_t IsPointInsidePlane(const Double_t *point, Double_t *vert, Double_t *norm) const;
   // Double_t SafetyToSector(const Double_t *point, Int_t iz, Double_t safmin, Bool_t in);
   // void SetIz(Int_t iz);
   // void SetSeg(Int_t iseg);

public:
   // constructors
   TGeoVGXtru();
   TGeoVGXtru(Int_t nz);
   TGeoVGXtru(Double_t *param);
   // destructor
   ~TGeoVGXtru() override;

   // methods
   // Function derived TGeoShape/TGeoBBox not present in TGeoVGAdapter
   // (not relevant to navigation)
   void ComputeBBox() override;
   Bool_t DefinePolygon(Int_t nvert, const Double_t *xv, const Double_t *yv);
   virtual void DefineSection(Int_t snum, Double_t z, Double_t x0 = 0., Double_t y0 = 0., Double_t scale = 1.);
   // Int_t DistancetoPrimitive(Int_t px, Int_t py) override;
   void DrawPolygon(Option_t *option = "");
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   //   virtual Int_t         GetByteCount() const {return 60+12*fNz;}

   Int_t GetNz() const { return fNz; }
   Int_t GetNvert() const { return fNvert; }
   Double_t GetX(Int_t i) const { return (i < fNvert && i > -1 && fX) ? fX[i] : -1.0E10; }
   Double_t GetY(Int_t i) const { return (i < fNvert && i > -1 && fY) ? fY[i] : -1.0E10; }
   Double_t GetXOffset(Int_t i) const { return (i < fNz && i > -1 && fX0) ? fX0[i] : 0.0; }
   Double_t GetYOffset(Int_t i) const { return (i < fNz && i > -1 && fY0) ? fY0[i] : 0.0; }
   Double_t GetScale(Int_t i) const { return (i < fNz && i > -1 && fScale) ? fScale[i] : 1.0; }
   Double_t *GetZ() const { return fZ; }
   Double_t GetZ(Int_t ipl) const;
   Double_t &Z(Int_t ipl) { return fZ[ipl]; }

   TGeoShape *GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const override { return nullptr; }
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNmeshVertices() const override;
   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetCurrentZ(Double_t z, Int_t iz);
   void SetCurrentVertices(Double_t x0, Double_t y0, Double_t scale);
   void SetDimensions(Double_t *param) override;
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // ClassDefOverride(TGeoVGXtru, 3) // extruded polygon class
};

#endif

#endif

