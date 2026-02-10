// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGPgon
#define ROOT_TGeoVGPgon

// TGeoVGPgon
//
// Class description:
//
// Wrapper class for TGeoPgon to make use of VecGeom Polyhedron.

#include "TGeoVGAdapter.h"

#include <iostream>

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#ifdef __ROOTCLING__
namespace vecgeom {
  // Providing a mock definition (Complete Type) for rootcling
  class UnplacedPolyhedron : public VUnplacedVolume {
    // nothing to be defined    
  };
}
#else
// The real build uses the actual VecGeom headers
#include <VecGeom/volumes/UnplacedPolyhedron.h>
#endif

class TGeoVGPgon final : public TGeoVGAdapter<vecgeom::UnplacedPolyhedron> {
  using Shape_t = vecgeom::UnplacedPolyhedron;
  using Base_t  = TGeoVGAdapter<vecgeom::UnplacedPolyhedron>;
private:
    // Static helper to provide "safe" values for the base class constructor
    static double* dummy() {
        // A minimal, valid polycone profile: two planes at z=0 and z=1
        // rmin = 0, rmax = 1. This avoids any "rmax < rmin" or "nz < 2" assertions.
        static double values[] = {0.0, 1.0}; 
        return values;
    }
public:
   struct ThreadData_t {
      Int_t *fIntBuffer;    //![fNedges+4] temporary int buffer array
      Double_t *fDblBuffer; //![fNedges+4] temporary double buffer array

      ThreadData_t();
      ~ThreadData_t();
   };
   ThreadData_t &GetThreadData() const;
   void ClearThreadData() const override;
   void CreateThreadData(Int_t nthreads) override;

protected:
   // data members
    // Data for construction
   Int_t fNz = 0;             // number of z planes (at least two)
   Double_t fPhi1 = 0;        // lower phi limit (converted to [0,2*pi)
   Double_t fDphi = 0;        // phi range
   Double_t *fRmin = nullptr; //[fNz] pointer to array of inner radii
   Double_t *fRmax = nullptr; //[fNz] pointer to array of outer radii
   Double_t *fZ = nullptr;    //[fNz] pointer to array of Z planes positions
   Bool_t fFullPhi = false;   //! Full phi range flag
   Double_t fC1 = 0;          //! Cosine of phi1
   Double_t fS1 = 0;          //! Sine of phi1
   Double_t fC2 = 0;          //! Cosine of phi1+dphi
   Double_t fS2 = 0;          //! Sine of phi1+dphi
   Double_t fCm = 0;          //! Cosine of (phi1+phi2)/2
   Double_t fSm = 0;          //! Sine of (phi1+phi2)/2
   Double_t fCdphi = 0;       //! Cosine of dphi

   // Pgon specific data
   Int_t fNedges;                                   // number of edges (at least one)
   mutable std::vector<ThreadData_t *> fThreadData; //! Navigation data per thread
   mutable Int_t fThreadSize;                       //! Size for the navigation data array
   mutable std::mutex fMutex;                       //! Mutex for thread data

       // The "actual" implementation
    vecgeom::UnplacedPolyhedron* fRealVolume = nullptr;

   // from TGeoPcon
   Bool_t HasInsideSurface() const;

   // internal utility methods
   // Int_t GetPhiCrossList(const Double_t *point, const Double_t *dir, Int_t istart, Double_t *sphi, Int_t *iphi,
   //                       Double_t stepmax = TGeoShape::Big()) const;
   // Bool_t IsCrossingSlice(const Double_t *point, const Double_t *dir, Int_t iphi, Double_t sstart, Int_t &ipl,
   //                        Double_t &snext, Double_t stepmax) const;
   // void LocatePhi(const Double_t *point, Int_t &ipsec) const;
   // Double_t Rpg(Double_t z, Int_t ipl, Bool_t inner, Double_t &a, Double_t &b) const;
   // Double_t Rproj(Double_t z, const Double_t *point, const Double_t *dir, Double_t cphi, Double_t sphi, Double_t &a,
   //                Double_t &b) const;
   // Bool_t SliceCrossing(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi, Double_t *sphi,
   //                      Double_t &snext, Double_t stepmax) const;
   // Bool_t SliceCrossingIn(const Double_t *point, const Double_t *dir, Int_t ipl, Int_t nphi, Int_t *iphi,
   //                        Double_t *sphi, Double_t &snext, Double_t stepmax) const;
   // Bool_t SliceCrossingZ(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi, Double_t *sphi,
   //                       Double_t &snext, Double_t stepmax) const;
   // Bool_t SliceCrossingInZ(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi, Double_t *sphi,
   //                         Double_t &snext, Double_t stepmax) const;
   void SetSegsAndPolsNoInside(TBuffer3D &buff) const;

   TGeoVGPgon(const TGeoVGPgon &) = delete;
   TGeoVGPgon &operator=(const TGeoVGPgon &) = delete;

public:
   // constructors
   TGeoVGPgon();
   TGeoVGPgon(Double_t phi, Double_t dphi, Int_t nedges, Int_t nz);
   TGeoVGPgon(const char *name, Double_t phi, Double_t dphi, Int_t nedges, Int_t nz);
   TGeoVGPgon(Double_t *params);
   // destructor
   ~TGeoVGPgon() override;

   // methods
   // methods
   // Method for construction
   virtual void DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax);

   // Methods from TGeoBBox
   void ComputeBBox() override;
   TGeoVolume *
   Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step) override;
   void GetBoundingCylinder(Double_t *param) const override;
   const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override;
   Int_t GetByteCount() const override { return 64 + 12 * fNz; }
   TGeoShape *GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const override { return nullptr; }
   void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override;
   Int_t GetNedges() const { return fNedges; }
   Int_t GetNmeshVertices() const override;
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override
   {
      return TGeoBBox::GetPointsOnSegments(npoints, array);
   }
   void InspectShape() const override;
   TBuffer3D *MakeBuffer3D() const override;
   Double_t SafetyToSegment(const Double_t *point, Int_t ipl, Int_t iphi, Bool_t in, Double_t safphi,
                            Double_t safmin = TGeoShape::Big()) const;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetDimensions(Double_t *param) override;
   void SetNedges(Int_t ne)
   {
      if (ne > 2)
         fNedges = ne;
   }
   void SetPoints(Double_t *points) const override;
   void SetPoints(Float_t *points) const override;
   void SetSegsAndPols(TBuffer3D &buff) const override;
   void Sizeof3D() const override;

   // Methods fron TGeoPcon
   Int_t GetNsegments() const { return fNedges; }
   Double_t GetPhi1() const { return fPhi1; }
   Double_t GetDphi() const { return fDphi; }
   Int_t GetNz() const { return fNz; }
   Double_t *GetRmin() const { return fRmin; }
   Double_t GetRmin(Int_t ipl) const;
   Double_t *GetRmax() const { return fRmax; }
   Double_t GetRmax(Int_t ipl) const;
   Double_t *GetZ() const { return fZ; }
   Double_t GetZ(Int_t ipl) const;

   Double_t &Phi1() { return fPhi1; }
   Double_t &Dphi() { return fDphi; }
   Double_t &Rmin(Int_t ipl) { return fRmin[ipl]; }
   Double_t &Rmax(Int_t ipl) { return fRmax[ipl]; }
   Double_t &Z(Int_t ipl) { return fZ[ipl]; }


   // Functions from GeoVGAdapter<UnplacedPolycone>
   // delegated to real volume
   //
   void ComputeNormal(const Double_t *point, const Double_t */*dir*/, Double_t *norm) const override {
      vecgeom::cxx::Vector3D<Double_t> vnorm;
      fRealVolume->Normal(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]), vnorm);
      norm[0] = vnorm.x();
      norm[1] = vnorm.y(), norm[2] = vnorm.z();
   }

   Bool_t Contains(const Double_t *point) const override {
      // return fRealVolume ? fRealVolume->Contains(p) : false;
      return fRealVolume->Contains(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]));
   }

   Double_t DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact = 1, Double_t step = TGeoShape::Big(),
                           Double_t *safe = nullptr) const override {
      Double_t dist = fRealVolume->DistanceToOut(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]),
                                              vecgeom::cxx::Vector3D<Double_t>(dir[0], dir[1], dir[2]), step);
      return ((dist < 0.) ? 0. : dist);
   }

   Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                            Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override {
      Double_t dist = fRealVolume->DistanceToIn(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]),
                                           vecgeom::cxx::Vector3D<Double_t>(dir[0], dir[1], dir[2]), step);
      return ((dist < 0.) ? 0. : dist);
   }

   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override {
      Double_t safety = (in) ? fRealVolume->SafetyToOut(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]))
                             : fRealVolume->SafetyToIn(vecgeom::cxx::Vector3D<Double_t>(point[0], point[1], point[2]));
      return ((safety < 0.) ? 0. : safety);
   }

    // VecGeom overridden methods ---------------------------------------------
   vecgeom::Precision
   DistanceToOut(U3Vector const& position, U3Vector const& direction,
                 vecgeom::Precision stepMax = vecgeom::kInfinityDbl) const override {
     std::cout << "Calling TGeoVGPgon DistanceToOut" << std::endl;
     return fRealVolume->DistanceToOut(position, direction, stepMax);
   }
   vecgeom::EnumInside
   Inside(U3Vector const& aPoint) const override {
     return fRealVolume->Inside(aPoint);
   }
   vecgeom::Precision
   DistanceToIn(U3Vector const& position, U3Vector const& direction,
                const vecgeom::Precision step_max = vecgeom::kInfinityDbl) const override {
     return fRealVolume->DistanceToIn(position, direction, step_max);
   }
   Bool_t Normal(U3Vector const& aPoint, U3Vector& aNormal) const override {
     return fRealVolume->Normal(aPoint, aNormal);
   }
   void Extent(U3Vector& aMin, U3Vector& aMax) const override {
     return fRealVolume->Extent(aMin, aMax);
   }
   U3Vector SamplePointOnSurface() const override {
     return fRealVolume->SamplePointOnSurface();
   }   

   // ClassDefOverride(TGeoVGPgon, 1) // polygone class
};

#endif // ROOT_USE_VECGEOM_SOLIDS

#endif
