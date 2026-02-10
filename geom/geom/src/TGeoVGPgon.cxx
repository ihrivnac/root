// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoVGPgon::Contains() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TGeoPgon.h"
#include "TGeoVGPgon.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoTube.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>


// ClassImp(TGeoVGPgon);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoVGPgon::ThreadData_t::ThreadData_t() : fIntBuffer(nullptr), fDblBuffer(nullptr) {}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVGPgon::ThreadData_t::~ThreadData_t()
{
   delete[] fIntBuffer;
   delete[] fDblBuffer;
}

////////////////////////////////////////////////////////////////////////////////

TGeoVGPgon::ThreadData_t &TGeoVGPgon::GetThreadData() const
{
   Int_t tid = TGeoManager::ThreadId();
   return *fThreadData[tid];
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVGPgon::ClearThreadData() const
{
   std::lock_guard<std::mutex> guard(fMutex);
   std::vector<ThreadData_t *>::iterator i = fThreadData.begin();
   while (i != fThreadData.end()) {
      delete *i;
      ++i;
   }
   fThreadData.clear();
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create thread data for n threads max.

void TGeoVGPgon::CreateThreadData(Int_t nthreads)
{
   if (fThreadSize)
      ClearThreadData();
   std::lock_guard<std::mutex> guard(fMutex);
   fThreadData.resize(nthreads);
   fThreadSize = nthreads;
   for (Int_t tid = 0; tid < nthreads; tid++) {
      if (fThreadData[tid] == nullptr) {
         fThreadData[tid] = new ThreadData_t;
         fThreadData[tid]->fIntBuffer = new Int_t[fNedges + 10];
         fThreadData[tid]->fDblBuffer = new Double_t[fNedges + 10];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoVGPgon::TGeoVGPgon()
   : Base_t("", 0., 360. * TMath::DegToRad(), 4, 2, dummy(), dummy(), dummy()),
     fNz(0),
     fPhi1(0.),
     fDphi(0.),
     fRmin(nullptr),
     fRmax(nullptr),
     fZ(nullptr),
     fFullPhi(kFALSE),
     fC1(0.),
     fS1(0.),
     fC2(0.),
     fS2(0.),
     fCm(0.),
     fSm(0.),
     fCdphi(0.),
     fNedges(0)
{
   SetShapeBit(TGeoShape::kGeoPgon);
   fNedges = 0;
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGPgon::TGeoVGPgon(Double_t phi, Double_t dphi, Int_t nedges, Int_t nz)
   : Base_t("", 0., 360. * TMath::DegToRad(), 4, 2, dummy(), dummy(), dummy()),
     fNz(nz),
     fPhi1(phi),
     fDphi(dphi),
     fRmin(nullptr),
     fRmax(nullptr),
     fZ(nullptr),
     fFullPhi(kFALSE),
     fC1(0.),
     fS1(0.),
     fC2(0.),
     fS2(0.),
     fCm(0.),
     fSm(0.),
     fCdphi(0.),
     fNedges(nedges)
{
   SetShapeBit(TGeoShape::kGeoPgon);
   // code from TGePcon
   while (fPhi1 < 0)
      fPhi1 += 360.;
   fRmin = new Double_t[nz];
   fRmax = new Double_t[nz];
   fZ = new Double_t[nz];
   memset(fRmin, 0, nz * sizeof(Double_t));
   memset(fRmax, 0, nz * sizeof(Double_t));
   memset(fZ, 0, nz * sizeof(Double_t));
   if (TGeoShape::IsSameWithinTolerance(fDphi, 360))
      fFullPhi = kTRUE;
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1 + fDphi;
   Double_t phim = 0.5 * (phi1 + phi2);
   fC1 = TMath::Cos(phi1 * TMath::DegToRad());
   fS1 = TMath::Sin(phi1 * TMath::DegToRad());
   fC2 = TMath::Cos(phi2 * TMath::DegToRad());
   fS2 = TMath::Sin(phi2 * TMath::DegToRad());
   fCm = TMath::Cos(phim * TMath::DegToRad());
   fSm = TMath::Sin(phim * TMath::DegToRad());
   fCdphi = TMath::Cos(0.5 * fDphi * TMath::DegToRad());
   // end of TGePcon code   
   fNedges = nedges;
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGPgon::TGeoVGPgon(const char *name, Double_t phi, Double_t dphi, Int_t nedges, Int_t nz)
   : Base_t(name, 0., 360. * TMath::DegToRad(), 4, 2, dummy(), dummy(), dummy()),
     fNz(nz),
     fPhi1(phi),
     fDphi(dphi),
     fRmin(nullptr),
     fRmax(nullptr),
     fZ(nullptr),
     fFullPhi(kFALSE),
     fC1(0.),
     fS1(0.),
     fC2(0.),
     fS2(0.),
     fCm(0.),
     fSm(0.),
     fCdphi(0.),
     fNedges(nedges)
{
   SetShapeBit(TGeoShape::kGeoPgon);
   // code from TGePcon
   while (fPhi1 < 0)
      fPhi1 += 360.;
   fRmin = new Double_t[nz];
   fRmax = new Double_t[nz];
   fZ = new Double_t[nz];
   memset(fRmin, 0, nz * sizeof(Double_t));
   memset(fRmax, 0, nz * sizeof(Double_t));
   memset(fZ, 0, nz * sizeof(Double_t));
   if (TGeoShape::IsSameWithinTolerance(fDphi, 360))
      fFullPhi = kTRUE;
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1 + fDphi;
   Double_t phim = 0.5 * (phi1 + phi2);
   fC1 = TMath::Cos(phi1 * TMath::DegToRad());
   fS1 = TMath::Sin(phi1 * TMath::DegToRad());
   fC2 = TMath::Cos(phi2 * TMath::DegToRad());
   fS2 = TMath::Sin(phi2 * TMath::DegToRad());
   fCm = TMath::Cos(phim * TMath::DegToRad());
   fSm = TMath::Sin(phim * TMath::DegToRad());
   fCdphi = TMath::Cos(0.5 * fDphi * TMath::DegToRad());
   // end of TGePcon code      
   fNedges = nedges;
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor in GEANT3 style
///  - param[0] = phi1
///  - param[1] = dphi
///  - param[2] = nedges
///  - param[3] = nz
///  - param[4] = z1
///  - param[5] = Rmin1
///  - param[6] = Rmax1
/// ...

TGeoVGPgon::TGeoVGPgon(Double_t *param)
   : Base_t("", 0., 360. * TMath::DegToRad(), 4, 2, dummy(), dummy(), dummy()),
     fNz(0),
     fPhi1(0.),
     fDphi(0.),
     fRmin(nullptr),
     fRmax(nullptr),
     fZ(nullptr),
     fFullPhi(kFALSE),
     fC1(0.),
     fS1(0.),
     fC2(0.),
     fS2(0.),
     fCm(0.),
     fSm(0.),
     fCdphi(0.),
     fNedges(0)
   {
   SetShapeBit(TGeoShape::kGeoPgon);
   SetDimensions(param);
   ComputeBBox();
   fThreadSize = 0;
   CreateThreadData(1);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGPgon::~TGeoVGPgon()
{
   ClearThreadData();
   if (fRmin) {
      delete[] fRmin;
      fRmin = nullptr;
   }
   if (fRmax) {
      delete[] fRmax;
      fRmax = nullptr;
   }
   if (fZ) {
      delete[] fZ;
      fZ = nullptr;
   }
   delete fRealVolume;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box for a polygone
/// Check if the sections are in increasing Z order

void TGeoVGPgon::ComputeBBox()
{
   for (Int_t isec = 0; isec < fNz - 1; isec++) {
      if (fZ[isec] > fZ[isec + 1]) {
         InspectShape();
         Fatal("ComputeBBox", "Wrong section order");
      }
   }
   // Check if the last sections are valid
   if (TMath::Abs(fZ[1] - fZ[0]) < TGeoShape::Tolerance() ||
       TMath::Abs(fZ[fNz - 1] - fZ[fNz - 2]) < TGeoShape::Tolerance()) {
      InspectShape();
      Fatal("ComputeBBox", "Shape %s at index %d: Not allowed first two or last two sections at same Z", GetName(),
            gGeoManager->GetListOfShapes()->IndexOf(this));
   }
   Double_t zmin = TMath::Min(fZ[0], fZ[fNz - 1]);
   Double_t zmax = TMath::Max(fZ[0], fZ[fNz - 1]);
   // find largest rmax an smallest rmin
   Double_t rmin, rmax;
   Double_t divphi = fDphi / fNedges;
   // find the radius of the outscribed circle
   rmin = fRmin[TMath::LocMin(fNz, fRmin)];
   rmax = fRmax[TMath::LocMax(fNz, fRmax)];
   rmax = rmax / TMath::Cos(0.5 * divphi * TMath::DegToRad());
   Double_t phi1 = fPhi1;
   Double_t phi2 = phi1 + fDphi;

   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax * TMath::Cos(phi1 * TMath::DegToRad());
   yc[0] = rmax * TMath::Sin(phi1 * TMath::DegToRad());
   xc[1] = rmax * TMath::Cos(phi2 * TMath::DegToRad());
   yc[1] = rmax * TMath::Sin(phi2 * TMath::DegToRad());
   xc[2] = rmin * TMath::Cos(phi1 * TMath::DegToRad());
   yc[2] = rmin * TMath::Sin(phi1 * TMath::DegToRad());
   xc[3] = rmin * TMath::Cos(phi2 * TMath::DegToRad());
   yc[3] = rmin * TMath::Sin(phi2 * TMath::DegToRad());

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];

   Double_t ddp = -phi1;
   if (ddp < 0)
      ddp += 360;
   if (ddp <= fDphi)
      xmax = rmax;
   ddp = 90 - phi1;
   if (ddp < 0)
      ddp += 360;
   if (ddp <= fDphi)
      ymax = rmax;
   ddp = 180 - phi1;
   if (ddp < 0)
      ddp += 360;
   if (ddp <= fDphi)
      xmin = -rmax;
   ddp = 270 - phi1;
   if (ddp < 0)
      ddp += 360;
   if (ddp <= fDphi)
      ymin = -rmax;
   fOrigin[0] = 0.5 * (xmax + xmin);
   fOrigin[1] = 0.5 * (ymax + ymin);
   fOrigin[2] = 0.5 * (zmax + zmin);
   fDX = 0.5 * (xmax - xmin);
   fDY = 0.5 * (ymax - ymin);
   fDZ = 0.5 * (zmax - zmin);
   SetShapeBit(kGeoClosedShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Defines z position of a section plane, rmin and rmax at this z. Sections
/// should be defined in increasing or decreasing Z order and the last section
/// HAS to be snum = fNz-1
/// From TGeoPcon

void TGeoVGPgon::DefineSection(Int_t snum, Double_t z, Double_t rmin, Double_t rmax)
{
   if ((snum < 0) || (snum >= fNz))
      return;
   fZ[snum] = z;
   fRmin[snum] = rmin;
   fRmax[snum] = rmax;
   if (rmin > rmax)
      Warning("DefineSection", "Shape %s: invalid rmin=%g rmax=%g", GetName(), rmin, rmax);
   if (snum == (fNz - 1)) {
      // Reorder sections in increasing Z order
      if (fZ[0] > fZ[snum]) {
         Int_t iz = 0;
         Int_t izi = fNz - 1;
         Double_t temp;
         while (iz < izi) {
            temp = fZ[iz];
            fZ[iz] = fZ[izi];
            fZ[izi] = temp;
            temp = fRmin[iz];
            fRmin[iz] = fRmin[izi];
            fRmin[izi] = temp;
            temp = fRmax[iz];
            fRmax[iz] = fRmax[izi];
            fRmax[izi] = temp;
            iz++;
            izi--;
         }
      }
      // Construct the VecGeom solid
      fRealVolume = new vecgeom::UnplacedPolyhedron(
        fPhi1 * TMath::DegToRad(), fDphi * TMath::DegToRad(), fNedges, fNz, fZ, fRmin, fRmax);
      ComputeBBox();
   }
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Locates index IPSEC of the phi sector containing POINT.

// void TGeoVGPgon::LocatePhi(const Double_t *point, Int_t &ipsec) const
// {
//    Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
//    while (phi < fPhi1)
//       phi += 360.;
//    ipsec = Int_t(fNedges * (phi - fPhi1) / fDphi); // [0, fNedges-1]
//    if (ipsec > fNedges - 1)
//       ipsec = -1; // in gap
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Returns lists of PGON phi crossings for a ray starting from POINT.

// Int_t TGeoVGPgon::GetPhiCrossList(const Double_t *point, const Double_t *dir, Int_t istart, Double_t *sphi, Int_t *iphi,
//                                 Double_t stepmax) const
// {
//    Double_t rxy, phi, cph, sph;
//    Int_t icrossed = 0;
//    if ((1. - TMath::Abs(dir[2])) < 1E-8) {
//       // ray is going parallel with Z
//       iphi[0] = istart;
//       sphi[0] = stepmax;
//       return 1;
//    }
//    Bool_t shootorig = (TMath::Abs(point[0] * dir[1] - point[1] * dir[0]) < 1E-8) ? kTRUE : kFALSE;
//    Double_t divphi = fDphi / fNedges;
//    if (shootorig) {
//       Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
//       if (rdotn > 0) {
//          sphi[0] = stepmax;
//          iphi[0] = istart;
//          return 1;
//       }
//       sphi[0] = TMath::Sqrt((point[0] * point[0] + point[1] * point[1]) / (1. - dir[2] * dir[2]));
//       iphi[0] = istart;
//       if (sphi[0] > stepmax) {
//          sphi[0] = stepmax;
//          return 1;
//       }
//       phi = TMath::ATan2(dir[1], dir[0]) * TMath::RadToDeg();
//       while (phi < fPhi1)
//          phi += 360.;
//       istart = Int_t((phi - fPhi1) / divphi);
//       if (istart > fNedges - 1)
//          istart = -1;
//       iphi[1] = istart;
//       sphi[1] = stepmax;
//       return 2;
//    }
//    Int_t incsec = Int_t(TMath::Sign(1., point[0] * dir[1] - point[1] * dir[0]));
//    Int_t ist;
//    if (istart < 0)
//       ist = (incsec > 0) ? 0 : fNedges;
//    else
//       ist = (incsec > 0) ? (istart + 1) : istart;
//    Bool_t crossing = kTRUE;
//    Bool_t gapdone = kFALSE;
//    divphi *= TMath::DegToRad();
//    Double_t phi1 = fPhi1 * TMath::DegToRad();
//    while (crossing) {
//       if (istart < 0)
//          gapdone = kTRUE;
//       phi = phi1 + ist * divphi;
//       cph = TMath::Cos(phi);
//       sph = TMath::Sin(phi);
//       crossing = IsCrossingSemiplane(point, dir, cph, sph, sphi[icrossed], rxy);
//       if (!crossing)
//          sphi[icrossed] = stepmax;
//       iphi[icrossed++] = istart;
//       if (crossing) {
//          if (sphi[icrossed - 1] > stepmax) {
//             sphi[icrossed - 1] = stepmax;
//             return icrossed;
//          }
//          if (istart < 0) {
//             istart = (incsec > 0) ? 0 : (fNedges - 1);
//          } else {
//             istart += incsec;
//             if (istart > fNedges - 1)
//                istart = (fDphi < 360.) ? (-1) : 0;
//             else if (istart < 0 && TGeoShape::IsSameWithinTolerance(fDphi, 360))
//                istart = fNedges - 1;
//          }
//          if (istart < 0) {
//             if (gapdone)
//                return icrossed;
//             ist = (incsec > 0) ? 0 : fNedges;
//          } else {
//             ist = (incsec > 0) ? (istart + 1) : istart;
//          }
//       }
//    }
//    return icrossed;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Performs ray propagation between Z segments.

// Bool_t TGeoVGPgon::SliceCrossingInZ(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi,
//                                   Double_t *stepphi, Double_t &snext, Double_t stepmax) const
// {
//    snext = 0.;
//    if (!nphi)
//       return kFALSE;
//    Int_t i;
//    Double_t rmin, rmax;
//    Double_t apg, bpg;
//    Double_t pt[3];
//    if (iphi[0] < 0 && nphi == 1)
//       return kFALSE;
//    // Get current Z segment
//    Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]);
//    if (ipl < 0 || ipl == fNz - 1)
//       return kFALSE;
//    if (TMath::Abs(point[2] - fZ[ipl]) < TGeoShape::Tolerance()) {
//       if (ipl < fNz - 2 && TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl + 1])) {
//          rmin = TMath::Min(fRmin[ipl], fRmin[ipl + 1]);
//          rmax = TMath::Max(fRmax[ipl], fRmax[ipl + 1]);
//       } else if (ipl > 1 && TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl - 1])) {
//          rmin = TMath::Min(fRmin[ipl], fRmin[ipl + 1]);
//          rmax = TMath::Max(fRmax[ipl], fRmax[ipl + 1]);
//       } else {
//          rmin = fRmin[ipl];
//          rmax = fRmax[ipl];
//       }
//    } else {
//       rmin = Rpg(point[2], ipl, kTRUE, apg, bpg);
//       rmax = Rpg(point[2], ipl, kFALSE, apg, bpg);
//    }
//    Int_t iphcrt;
//    Double_t divphi = TMath::DegToRad() * fDphi / fNedges;
//    Double_t rproj, ndot, dist;
//    Double_t phi1 = fPhi1 * TMath::DegToRad();
//    Double_t cosph, sinph;
//    Double_t snextphi = 0.;
//    Double_t step = 0;
//    Double_t phi;
//    memcpy(pt, point, 3 * sizeof(Double_t));
//    for (iphcrt = 0; iphcrt < nphi; iphcrt++) {
//       if (step > stepmax) {
//          snext = step;
//          return kFALSE;
//       }
//       if (iphi[iphcrt] < 0) {
//          snext = step;
//          return kTRUE;
//       }
//       // check crossing
//       snextphi = stepphi[iphcrt];
//       phi = phi1 + (iphi[iphcrt] + 0.5) * divphi;
//       cosph = TMath::Cos(phi);
//       sinph = TMath::Sin(phi);
//       rproj = pt[0] * cosph + pt[1] * sinph;
//       dist = TGeoShape::Big();
//       ndot = dir[0] * cosph + dir[1] * sinph;
//       if (!TGeoShape::IsSameWithinTolerance(ndot, 0)) {
//          dist = (ndot > 0) ? ((rmax - rproj) / ndot) : ((rmin - rproj) / ndot);
//          if (dist < 0)
//             dist = 0.;
//       }
//       if (dist < (snextphi - step)) {
//          snext = step + dist;
//          if (snext < stepmax)
//             return kTRUE;
//          return kFALSE;
//       }
//       step = snextphi;
//       for (i = 0; i < 3; i++)
//          pt[i] = point[i] + step * dir[i];
//    }
//    snext = step;
//    return kFALSE;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Performs ray propagation between Z segments.

// Bool_t TGeoVGPgon::SliceCrossingZ(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi, Double_t *stepphi,
//                                 Double_t &snext, Double_t stepmax) const
// {
//    if (!nphi)
//       return kFALSE;
//    Int_t i;
//    Double_t rmin, rmax;
//    Double_t apg, bpg;
//    Double_t pt[3];
//    if (iphi[0] < 0 && nphi == 1)
//       return kFALSE;
//    // Get current Z segment
//    Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]);
//    if (ipl < 0 || ipl == fNz - 1)
//       return kFALSE;
//    if (TMath::Abs(point[2] - fZ[ipl]) < TGeoShape::Tolerance()) {
//       if (ipl < fNz - 2 && TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl + 1])) {
//          rmin = TMath::Min(fRmin[ipl], fRmin[ipl + 1]);
//          rmax = TMath::Max(fRmax[ipl], fRmax[ipl + 1]);
//       } else if (ipl > 1 && TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl - 1])) {
//          rmin = TMath::Min(fRmin[ipl], fRmin[ipl + 1]);
//          rmax = TMath::Max(fRmax[ipl], fRmax[ipl + 1]);
//       } else {
//          rmin = fRmin[ipl];
//          rmax = fRmax[ipl];
//       }
//    } else {
//       rmin = Rpg(point[2], ipl, kTRUE, apg, bpg);
//       rmax = Rpg(point[2], ipl, kFALSE, apg, bpg);
//    }
//    Int_t iphcrt;
//    Double_t divphi = TMath::DegToRad() * fDphi / fNedges;
//    Double_t rproj, ndot, dist;
//    Double_t phi1 = fPhi1 * TMath::DegToRad();
//    Double_t cosph, sinph;
//    Double_t snextphi = 0.;
//    Double_t step = 0;
//    Double_t phi;
//    memcpy(pt, point, 3 * sizeof(Double_t));
//    for (iphcrt = 0; iphcrt < nphi; iphcrt++) {
//       if (step > stepmax)
//          return kFALSE;
//       snextphi = stepphi[iphcrt];
//       if (iphi[iphcrt] < 0) {
//          if (iphcrt == nphi - 1)
//             return kFALSE;
//          if (snextphi > stepmax)
//             return kFALSE;
//          for (i = 0; i < 3; i++)
//             pt[i] = point[i] + snextphi * dir[i];
//          phi = phi1 + (iphi[iphcrt + 1] + 0.5) * divphi;
//          cosph = TMath::Cos(phi);
//          sinph = TMath::Sin(phi);
//          rproj = pt[0] * cosph + pt[1] * sinph;
//          if (rproj < rmin || rproj > rmax) {
//             step = snextphi;
//             continue;
//          }
//          snext = snextphi;
//          return kTRUE;
//       }
//       // check crossing
//       phi = phi1 + (iphi[iphcrt] + 0.5) * divphi;
//       cosph = TMath::Cos(phi);
//       sinph = TMath::Sin(phi);
//       rproj = pt[0] * cosph + pt[1] * sinph;
//       dist = TGeoShape::Big();
//       ndot = dir[0] * cosph + dir[1] * sinph;
//       if (rproj < rmin) {
//          dist = (ndot > 0) ? ((rmin - rproj) / ndot) : TGeoShape::Big();
//       } else {
//          dist = (ndot < 0) ? ((rmax - rproj) / ndot) : TGeoShape::Big();
//       }
//       if (dist < 1E10) {
//          snext = step + dist;
//          if (snext < stepmax)
//             return kTRUE;
//       }
//       step = snextphi;
//       for (i = 0; i < 3; i++)
//          pt[i] = point[i] + step * dir[i];
//    }
//    return kFALSE;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Check boundary crossing inside phi slices. Return distance snext to first crossing
// /// if smaller than stepmax.
// /// Protection in case point is in phi gap or close to phi boundaries and exiting

// Bool_t TGeoVGPgon::SliceCrossingIn(const Double_t *point, const Double_t *dir, Int_t ipl, Int_t nphi, Int_t *iphi,
//                                  Double_t *stepphi, Double_t &snext, Double_t stepmax) const
// {
//    snext = 0.;
//    if (!nphi)
//       return kFALSE;
//    Int_t i;
//    Int_t iphstart = 0;
//    Double_t pt[3];
//    if (iphi[0] < 0) {
//       if (stepphi[0] > TGeoShape::Tolerance())
//          return kFALSE;
//       iphstart = 1;
//    }
//    if (nphi > 1 && iphi[1] < 0 && stepphi[0] < TGeoShape::Tolerance()) {
//       snext = stepphi[0];
//       return kTRUE;
//    }
//    // Get current Z segment
//    Double_t snextphi = 0.;
//    Double_t step = 0;
//    Int_t incseg = (dir[2] > 0) ? 1 : -1; // dir[2] is never 0 here
//    // Compute the projected radius from starting point
//    Int_t iplstart = ipl;
//    Int_t iphcrt = 0;
//    Double_t apr = TGeoShape::Big(), bpr = 0, db = 0;
//    Double_t rpg = 0, rnew = 0, znew = 0;
//    Double_t rpgin = 0, rpgout = 0, apgin = 0, apgout = 0, bpgin = 0, bpgout = 0;
//    Double_t divphi = TMath::DegToRad() * fDphi / fNedges;
//    Double_t phi1 = fPhi1 * TMath::DegToRad();
//    Double_t phi = 0, dz = 0;
//    Double_t cosph = 0, sinph = 0;
//    Double_t distz = 0, distr = 0, din = 0, dout = 0;
//    Double_t invdir = 1. / dir[2];
//    memcpy(pt, point, 3 * sizeof(Double_t));
//    for (iphcrt = iphstart; iphcrt < nphi; iphcrt++) {
//       // check if step to current checked slice is too big
//       if (step > stepmax) {
//          snext = step;
//          return kFALSE;
//       }
//       if (iphi[iphcrt] < 0) {
//          snext = snextphi;
//          return kTRUE;
//       }
//       snextphi = stepphi[iphcrt];
//       phi = phi1 + (iphi[iphcrt] + 0.5) * divphi;
//       cosph = TMath::Cos(phi);
//       sinph = TMath::Sin(phi);
//       Double_t rproj = Rproj(pt[2], pt, dir, cosph, sinph, apr, bpr);
//       // compute distance to next Z plane
//       while (ipl >= 0 && ipl < fNz - 1) {
//          din = dout = TGeoShape::Big();
//          // dist to last boundary of current segment according dir
//          distz = (fZ[ipl + ((1 + incseg) >> 1)] - pt[2]) * invdir;
//          // length of current segment
//          dz = fZ[ipl + 1] - fZ[ipl];
//          if (dz < TGeoShape::Tolerance()) {
//             rnew = apr + bpr * fZ[ipl];
//             rpg = (rnew - fRmin[ipl]) * (rnew - fRmin[ipl + 1]);
//             if (rpg <= 0)
//                din = distz;
//             rpg = (rnew - fRmax[ipl]) * (rnew - fRmax[ipl + 1]);
//             if (rpg <= 0)
//                dout = distz;
//             distr = TMath::Min(din, dout);
//          } else {
//             rpgin = Rpg(pt[2], ipl, kTRUE, apgin, bpgin);
//             db = bpgin - bpr;
//             if (TMath::Abs(db) > TGeoShape::Tolerance()) {
//                znew = (apr - apgin) / db;
//                din = (znew - pt[2]) * invdir;
//             }
//             rpgout = Rpg(pt[2], ipl, kFALSE, apgout, bpgout);
//             db = bpgout - bpr;
//             if (TMath::Abs(db) > TGeoShape::Tolerance()) {
//                znew = (apr - apgout) / db;
//                dout = (znew - pt[2]) * invdir;
//             }
//             // protection for the first segment
//             Double_t dinp = (din > TMath::Abs(snext - TGeoShape::Tolerance())) ? din : TGeoShape::Big();
//             Double_t doutp = (dout > TMath::Abs(snext - TGeoShape::Tolerance())) ? dout : TGeoShape::Big();
//             distr = TMath::Min(dinp, doutp);
//             if (iphcrt == iphstart && ipl == iplstart) {
//                if (rproj < rpgin + 1.E-8) {
//                   Double_t ndotd = dir[0] * cosph + dir[1] * sinph + dir[2] * (fRmin[ipl] - fRmin[ipl + 1]) / dz;
//                   if (ndotd < 0) {
//                      snext = (din < 0) ? step : (step + din);
//                      return kTRUE;
//                   } else {
//                      // Ignore din
//                      din = -TGeoShape::Big();
//                   }
//                   distr = TMath::Max(din, dout);
//                   if (distr < TGeoShape::Tolerance())
//                      distr = TGeoShape::Big();
//                } else if (rproj > rpgout - 1.E-8) {
//                   Double_t ndotd = dir[0] * cosph + dir[1] * sinph + dir[2] * (fRmax[ipl] - fRmax[ipl + 1]) / dz;
//                   if (ndotd > 0) {
//                      snext = (dout < 0) ? step : (step + dout);
//                      return kTRUE;
//                   } else {
//                      // Ignore dout
//                      dout = -TGeoShape::Big();
//                   }
//                   distr = TMath::Max(din, dout);
//                   if (distr < TGeoShape::Tolerance())
//                      distr = TGeoShape::Big();
//                }
//             }
//          }
//          if (distr < snext - TGeoShape::Tolerance())
//             distr = TGeoShape::Big();
//          if (snextphi < step + TMath::Min(distz, distr)) {
//             for (i = 0; i < 3; i++)
//                pt[i] = point[i] + snextphi * dir[i];
//             step = snextphi;
//             snext = 0.0;
//             break;
//          }
//          if (distr <= distz + TGeoShape::Tolerance()) {
//             step += distr;
//             snext = step;
//             return (step > stepmax) ? kFALSE : kTRUE;
//          }
//          // we have crossed a Z boundary
//          snext = distz;
//          if ((ipl + incseg < 0) || (ipl + incseg > fNz - 2)) {
//             // it was the last boundary
//             step += distz;
//             snext = step;
//             return (step > stepmax) ? kFALSE : kTRUE;
//          }
//          ipl += incseg;
//       } // end loop Z
//    }    // end loop phi
//    snext = TGeoShape::Big();
//    return kFALSE;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Check boundary crossing inside phi slices. Return distance snext to first crossing
// /// if smaller than stepmax.

// Bool_t TGeoVGPgon::SliceCrossing(const Double_t *point, const Double_t *dir, Int_t nphi, Int_t *iphi, Double_t *stepphi,
//                                Double_t &snext, Double_t stepmax) const
// {
//    if (!nphi)
//       return kFALSE;
//    Int_t i;
//    Double_t pt[3];
//    if (iphi[0] < 0 && nphi == 1)
//       return kFALSE;

//    Double_t snextphi = 0.;
//    Double_t step = 0;
//    // Get current Z segment
//    Int_t incseg = (dir[2] > 0) ? 1 : -1; // dir[2] is never 0 here
//    Int_t ipl = TMath::BinarySearch(fNz, fZ, point[2]);
//    if (ipl < 0) {
//       ipl = 0; // this should never happen
//       if (incseg < 0)
//          return kFALSE;
//    } else {
//       if (ipl == fNz - 1) {
//          ipl = fNz - 2; // nor this
//          if (incseg > 0)
//             return kFALSE;
//       } else {
//          if (TMath::Abs(point[2] - fZ[ipl]) < TGeoShape::Tolerance()) {
//             // we are at the sector edge, but never inside the pgon
//             if ((ipl + incseg) < 0 || (ipl + incseg) > fNz - 1)
//                return kFALSE;
//             if (TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl + incseg]))
//                ipl += incseg;
//             // move to next clean segment if downwards
//             if (incseg < 0) {
//                if (TGeoShape::IsSameWithinTolerance(fZ[ipl], fZ[ipl + 1]))
//                   ipl--;
//             }
//          }
//       }
//    }
//    // Compute the projected radius from starting point
//    Int_t iphcrt;
//    Double_t apg, bpg;
//    Double_t rpgin;
//    Double_t rpgout;
//    Double_t divphi = TMath::DegToRad() * fDphi / fNedges;
//    Double_t phi1 = fPhi1 * TMath::DegToRad();
//    Double_t phi;
//    Double_t cosph, sinph;
//    Double_t rproj;
//    memcpy(pt, point, 3 * sizeof(Double_t));
//    for (iphcrt = 0; iphcrt < nphi; iphcrt++) {
//       // check if step to current checked slice is too big
//       if (step > stepmax)
//          return kFALSE;
//       // jump over the dead sector
//       snextphi = stepphi[iphcrt];
//       if (iphi[iphcrt] < 0) {
//          if (iphcrt == nphi - 1)
//             return kFALSE;
//          if (snextphi > stepmax)
//             return kFALSE;
//          for (i = 0; i < 3; i++)
//             pt[i] = point[i] + snextphi * dir[i];
//          // we have a new z, so check again iz
//          if (incseg > 0) {
//             // loop z planes
//             while (pt[2] > fZ[ipl + 1]) {
//                ipl++;
//                if (ipl > fNz - 2)
//                   return kFALSE;
//             }
//          } else {
//             while (pt[2] < fZ[ipl]) {
//                ipl--;
//                if (ipl < 0)
//                   return kFALSE;
//             }
//          }
//          // check if we have a crossing when entering new sector
//          rpgin = Rpg(pt[2], ipl, kTRUE, apg, bpg);
//          rpgout = Rpg(pt[2], ipl, kFALSE, apg, bpg);
//          phi = phi1 + (iphi[iphcrt + 1] + 0.5) * divphi;
//          cosph = TMath::Cos(phi);
//          sinph = TMath::Sin(phi);

//          rproj = pt[0] * cosph + pt[1] * sinph;
//          if (rproj < rpgin || rproj > rpgout) {
//             step = snextphi;
//             continue;
//          }
//          snext = snextphi;
//          return kTRUE;
//       }
//       if (IsCrossingSlice(point, dir, iphi[iphcrt], step, ipl, snext, TMath::Min(snextphi, stepmax)))
//          return kTRUE;
//       step = snextphi;
//    }
//    return kFALSE;
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Check crossing of a given pgon slice, from a starting point inside the slice

// Bool_t TGeoVGPgon::IsCrossingSlice(const Double_t *point, const Double_t *dir, Int_t iphi, Double_t sstart, Int_t &ipl,
//                                  Double_t &snext, Double_t stepmax) const
// {
//    if (ipl < 0 || ipl > fNz - 2)
//       return kFALSE;
//    if (sstart > stepmax)
//       return kFALSE;
//    Double_t pt[3];
//    memcpy(pt, point, 3 * sizeof(Double_t));
//    if (sstart > 0)
//       for (Int_t i = 0; i < 3; i++)
//          pt[i] += sstart * dir[i];
//    stepmax -= sstart;
//    Double_t step;
//    Int_t incseg = (dir[2] > 0) ? 1 : -1;
//    Double_t invdir = 1. / dir[2];
//    Double_t divphi = TMath::DegToRad() * fDphi / fNedges;
//    Double_t phi = fPhi1 * TMath::DegToRad() + (iphi + 0.5) * divphi;
//    Double_t cphi = TMath::Cos(phi);
//    Double_t sphi = TMath::Sin(phi);
//    Double_t apr = TGeoShape::Big();
//    Double_t bpr = 0.;
//    Rproj(pt[2], point, dir, cphi, sphi, apr, bpr);
//    Double_t dz;
//    // loop segments
//    Int_t icrtseg = ipl;
//    Int_t isegstart = ipl;
//    Int_t iseglast = (incseg > 0) ? (fNz - 1) : -1;
//    Double_t din, dout, rdot, rnew, apg, bpg, db, znew;

//    for (ipl = isegstart; ipl != iseglast; ipl += incseg) {
//       step = (fZ[ipl + 1 - ((1 + incseg) >> 1)] - pt[2]) * invdir;
//       if (step > 0) {
//          if (step > stepmax) {
//             ipl = icrtseg;
//             return kFALSE;
//          }
//          icrtseg = ipl;
//       }
//       din = dout = TGeoShape::Big();
//       dz = fZ[ipl + 1] - fZ[ipl];

//       //      rdot = (rproj-fRmin[ipl])*dz - (pt[2]-fZ[ipl])*(fRmin[ipl+1]-fRmin[ipl]);
//       if (TGeoShape::IsSameWithinTolerance(dz, 0))
//          rdot = dir[2] * TMath::Sign(1., fRmin[ipl] - fRmin[ipl + 1]);
//       else
//          rdot = dir[0] * cphi + dir[1] * sphi + dir[2] * (fRmin[ipl] - fRmin[ipl + 1]) / dz;
//       if (rdot > 0) {
//          // inner surface visible ->check crossing
//          //         printf("   inner visible\n");
//          if (TGeoShape::IsSameWithinTolerance(dz, 0)) {
//             rnew = apr + bpr * fZ[ipl];
//             Double_t rpg = (rnew - fRmin[ipl]) * (rnew - fRmin[ipl + 1]);
//             if (rpg <= 0)
//                din = (fZ[ipl] - pt[2]) * invdir;
//          } else {
//             Rpg(pt[2], ipl, kTRUE, apg, bpg);
//             db = bpg - bpr;
//             if (!TGeoShape::IsSameWithinTolerance(db, 0)) {
//                znew = (apr - apg) / db;
//                if (znew > fZ[ipl] && znew < fZ[ipl + 1]) {
//                   din = (znew - pt[2]) * invdir;
//                   if (din < 0)
//                      din = TGeoShape::Big();
//                }
//             }
//          }
//       }
//       //      printf("   din=%f\n", din);
//       //      rdot = (rproj-fRmax[ipl])*dz - (pt[2]-fZ[ipl])*(fRmax[ipl+1]-fRmax[ipl]);
//       if (TGeoShape::IsSameWithinTolerance(dz, 0))
//          rdot = dir[2] * TMath::Sign(1., fRmax[ipl] - fRmax[ipl + 1]);
//       else
//          rdot = dir[0] * cphi + dir[1] * sphi + dir[2] * (fRmax[ipl] - fRmax[ipl + 1]) / dz;
//       if (rdot < 0) {
//          //         printf("   outer visible\n");
//          // outer surface visible ->check crossing
//          if (TGeoShape::IsSameWithinTolerance(dz, 0)) {
//             rnew = apr + bpr * fZ[ipl];
//             Double_t rpg = (rnew - fRmax[ipl]) * (rnew - fRmax[ipl + 1]);
//             if (rpg <= 0)
//                dout = (fZ[ipl] - pt[2]) * invdir;
//          } else {
//             Rpg(pt[2], ipl, kFALSE, apg, bpg);
//             db = bpg - bpr;
//             if (!TGeoShape::IsSameWithinTolerance(db, 0)) {
//                znew = (apr - apg) / db;
//                if (znew > fZ[ipl] && znew < fZ[ipl + 1])
//                   dout = (znew - pt[2]) * invdir;
//                if (dout < 0)
//                   dout = TGeoShape::Big();
//             }
//          }
//       }
//       //      printf("   dout=%f\n", dout);
//       step = TMath::Min(din, dout);
//       if (step < 1E10) {
//          // there is a crossing within this segment
//          if (step > stepmax) {
//             ipl = icrtseg;
//             return kFALSE;
//          }
//          snext = sstart + step;
//          return kTRUE;
//       }
//    }
//    ipl = icrtseg;
//    return kFALSE;
// }

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGPgon::GetBoundingCylinder(Double_t *param) const
{
   param[0] = fRmin[0]; // Rmin
   param[1] = fRmax[0]; // Rmax
   for (Int_t i = 1; i < fNz; i++) {
      if (fRmin[i] < param[0])
         param[0] = fRmin[i];
      if (fRmax[i] > param[1])
         param[1] = fRmax[i];
   }
   Double_t divphi = fDphi / fNedges;
   param[1] /= TMath::Cos(0.5 * divphi * TMath::DegToRad());
   param[0] *= param[0];
   param[1] *= param[1];
   if (TGeoShape::IsSameWithinTolerance(fDphi, 360)) {
      param[2] = 0.;
      param[3] = 360.;
      return;
   }
   param[2] = (fPhi1 < 0) ? (fPhi1 + 360.) : fPhi1; // Phi1
   param[3] = param[2] + fDphi;                     // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this polygone shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. Phi divisions are
/// allowed only if nedges%ndiv=0 and create polygone "segments" with nedges/ndiv edges.
/// Z divisions can be performed if the divided range is in between two consecutive Z planes.
/// In case a wrong division axis is supplied, returns pointer to volume that was divided.

TGeoVolume *
TGeoVGPgon::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   //   printf("Dividing %s : nz=%d nedges=%d phi1=%g dphi=%g (ndiv=%d iaxis=%d start=%g step=%g)\n",
   //          voldiv->GetName(), fNz, fNedges, fPhi1, fDphi, ndiv, iaxis, start, step);
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Int_t nedges = fNedges;
   Double_t zmin = start;
   Double_t zmax = start + ndiv * step;
   Int_t isect = -1;
   Int_t is, id, ipl;
   switch (iaxis) {
   case 1: //---                R division
      Error("Divide", "makes no sense dividing a pgon on radius");
      return nullptr;
   case 2: //---                Phi division
      if (fNedges % ndiv) {
         Error("Divide", "ndiv should divide number of pgon edges");
         return nullptr;
      }
      nedges = fNedges / ndiv;
      finder = new TGeoPatternCylPhi(voldiv, ndiv, start, start + ndiv * step);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGPgon(-step / 2, step, nedges, fNz);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      for (is = 0; is < fNz; is++)
         ((TGeoVGPgon *)shape)->DefineSection(is, fZ[is], fRmin[is], fRmax[is]);
      opt = "Phi";
      for (id = 0; id < ndiv; id++) {
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 3: // ---                Z division
      // find start plane
      for (ipl = 0; ipl < fNz - 1; ipl++) {
         if (start < fZ[ipl])
            continue;
         else {
            if ((start + ndiv * step) > fZ[ipl + 1])
               continue;
         }
         isect = ipl;
         zmin = fZ[isect];
         zmax = fZ[isect + 1];
         break;
      }
      if (isect < 0) {
         Error("Divide", "cannot divide pcon on Z if divided region is not between 2 consecutive planes");
         return nullptr;
      }
      finder = new TGeoPatternZ(voldiv, ndiv, start, start + ndiv * step);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      opt = "Z";
      for (id = 0; id < ndiv; id++) {
         Double_t z1 = start + id * step;
         Double_t z2 = start + (id + 1) * step;
         Double_t rmin1 = (fRmin[isect] * (zmax - z1) - fRmin[isect + 1] * (zmin - z1)) / (zmax - zmin);
         Double_t rmax1 = (fRmax[isect] * (zmax - z1) - fRmax[isect + 1] * (zmin - z1)) / (zmax - zmin);
         Double_t rmin2 = (fRmin[isect] * (zmax - z2) - fRmin[isect + 1] * (zmin - z2)) / (zmax - zmin);
         Double_t rmax2 = (fRmax[isect] * (zmax - z2) - fRmax[isect + 1] * (zmin - z2)) / (zmax - zmin);
         shape = new TGeoVGPgon(fPhi1, fDphi, nedges, 2);
         ((TGeoVGPgon *)shape)->DefineSection(0, -step / 2, rmin1, rmax1);
         ((TGeoVGPgon *)shape)->DefineSection(1, step / 2, rmin2, rmax2);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   default: Error("Divide", "Wrong axis type for division"); return nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Rmin for Z segment IPL.
/// From TGeoPcon

Double_t TGeoVGPgon::GetRmin(Int_t ipl) const
{
   if (ipl < 0 || ipl > (fNz - 1)) {
      Error("GetRmin", "ipl=%i out of range (0,%i) in shape %s", ipl, fNz - 1, GetName());
      return 0.;
   }
   return fRmin[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Rmax for Z segment IPL.
/// From TGeoPcon

Double_t TGeoVGPgon::GetRmax(Int_t ipl) const
{
   if (ipl < 0 || ipl > (fNz - 1)) {
      Error("GetRmax", "ipl=%i out of range (0,%i) in shape %s", ipl, fNz - 1, GetName());
      return 0.;
   }
   return fRmax[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Returns Z for segment IPL.
/// From TGeoPcon

Double_t TGeoVGPgon::GetZ(Int_t ipl) const
{
   if (ipl < 0 || ipl > (fNz - 1)) {
      Error("GetZ", "ipl=%i out of range (0,%i) in shape %s", ipl, fNz - 1, GetName());
      return 0.;
   }
   return fZ[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Inspect the PGON parameters.

void TGeoVGPgon::InspectShape() const
{
   printf("*** Shape %s: TGeoVGPgon ***\n", GetName());
   printf("    Nedges = %i\n", fNedges);
   printf("    Nz    = %i\n", fNz);
   printf("    phi1  = %11.5f\n", fPhi1);
   printf("    dphi  = %11.5f\n", fDphi);
   for (Int_t ipl = 0; ipl < fNz; ipl++)
      printf("     plane %i: z=%11.5f Rmin=%11.5f Rmax=%11.5f\n", ipl, fZ[ipl], fRmin[ipl], fRmax[ipl]);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGPgon::MakeBuffer3D() const
{
   Int_t nbPnts, nbSegs, nbPols;
   GetMeshNumbers(nbPnts, nbSegs, nbPols);

   if (nbPnts <= 0)
      return nullptr;

   TBuffer3D *buff =
      new TBuffer3D(TBuffer3DTypes::kGeneric, nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols);
   if (buff) {
      SetPoints(buff->fPnts);
      SetSegsAndPols(*buff);
   }

   return buff;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons.

void TGeoVGPgon::SetSegsAndPols(TBuffer3D &buff) const
{
   if (!HasInsideSurface()) {
      SetSegsAndPolsNoInside(buff);
      return;
   }

   Int_t i, j;
   const Int_t n = GetNedges() + 1;
   Int_t nz = GetNz();
   if (nz < 2)
      return;
   Int_t nbPnts = nz * 2 * n;
   if (nbPnts <= 0)
      return;
   Double_t dphi = GetDphi();
   Bool_t specialCase = TGeoShape::IsSameWithinTolerance(dphi, 360);

   Int_t c = GetBasicColor();

   Int_t indx = 0, indx2, k;

   // inside & outside circles, number of segments: 2*nz*(n-1)
   //             special case number of segments: 2*nz*n
   for (i = 0; i < nz * 2; i++) {
      indx2 = i * n;
      for (j = 1; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j - 1;
         buff.fSegs[indx++] = indx2 + j;
      }
      if (specialCase) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j - 1;
         buff.fSegs[indx++] = indx2;
      }
   }

   // bottom & top lines, number of segments: 2*n
   for (i = 0; i < 2; i++) {
      indx2 = i * (nz - 1) * 2 * n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j;
         buff.fSegs[indx++] = indx2 + n + j;
      }
   }

   // inside & outside cylinders, number of segments: 2*(nz-1)*n
   for (i = 0; i < (nz - 1); i++) {
      // inside cylinder
      indx2 = i * n * 2;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = indx2 + j;
         buff.fSegs[indx++] = indx2 + n * 2 + j;
      }
      // outside cylinder
      indx2 = i * n * 2 + n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c + 3;
         buff.fSegs[indx++] = indx2 + j;
         buff.fSegs[indx++] = indx2 + n * 2 + j;
      }
   }

   // left & right sections, number of segments: 2*(nz-2)
   //          special case number of segments: 0
   if (!specialCase) {
      for (i = 1; i < (nz - 1); i++) {
         for (j = 0; j < 2; j++) {
            buff.fSegs[indx++] = c;
            buff.fSegs[indx++] = 2 * i * n + j * (n - 1);
            buff.fSegs[indx++] = (2 * i + 1) * n + j * (n - 1);
         }
      }
   }

   Int_t m = n - 1 + (specialCase ? 1 : 0);
   indx = 0;

   // bottom & top, number of polygons: 2*(n-1)
   // special case number of polygons: 2*n
   i = 0;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c + 3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2 * nz * m + i * n + j;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n + j + 1;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + j;
   }
   if (specialCase) {
      buff.fPols[indx++] = c + 3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = 2 * nz * m + i * n + j;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + j;
   }
   i = 1;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c + 3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n + j + 1;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n + j;
   }
   if (specialCase) {
      buff.fPols[indx++] = c + 3;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n;
      buff.fPols[indx++] = i * (nz * 2 - 2) * m + m + j;
      buff.fPols[indx++] = 2 * nz * m + i * n + j;
   }

   // inside & outside, number of polygons: (nz-1)*2*(n-1)
   for (k = 0; k < (nz - 1); k++) {
      i = 0;
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c + i;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j + 1;
         buff.fPols[indx++] = (2 * k + i * 1 + 2) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j;
         buff.fPols[indx++] = (2 * k + i * 1) * m + j;
      }
      if (specialCase) {
         buff.fPols[indx++] = c + i;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n;
         buff.fPols[indx++] = (2 * k + i * 1 + 2) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j;
         buff.fPols[indx++] = (2 * k + i * 1) * m + j;
      }
      i = 1;
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c + i;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = (2 * k + i * 1) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j;
         buff.fPols[indx++] = (2 * k + i * 1 + 2) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j + 1;
      }
      if (specialCase) {
         buff.fPols[indx++] = c + i;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = (2 * k + i * 1) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n + j;
         buff.fPols[indx++] = (2 * k + i * 1 + 2) * m + j;
         buff.fPols[indx++] = nz * 2 * m + (2 * k + i * 1 + 2) * n;
      }
   }

   // left & right sections, number of polygons: 2*(nz-1)
   //          special case number of polygons: 0
   if (!specialCase) {
      indx2 = nz * 2 * (n - 1);
      for (k = 0; k < (nz - 1); k++) {
         buff.fPols[indx++] = c + 2;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = k == 0 ? indx2 : indx2 + 2 * nz * n + 2 * (k - 1);
         buff.fPols[indx++] = indx2 + 2 * (k + 1) * n;
         buff.fPols[indx++] = indx2 + 2 * nz * n + 2 * k;
         buff.fPols[indx++] = indx2 + (2 * k + 3) * n;

         buff.fPols[indx++] = c + 2;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = k == 0 ? indx2 + n - 1 : indx2 + 2 * nz * n + 2 * (k - 1) + 1; // a
         buff.fPols[indx++] = indx2 + (2 * k + 3) * n + n - 1;                               // d
         buff.fPols[indx++] = indx2 + 2 * nz * n + 2 * k + 1;                                // c
         buff.fPols[indx++] = indx2 + 2 * (k + 1) * n + n - 1;                               // b
      }
      buff.fPols[indx - 8] = indx2 + n;
      buff.fPols[indx - 2] = indx2 + 2 * n - 1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons, when no inner surface exists

void TGeoVGPgon::SetSegsAndPolsNoInside(TBuffer3D &buff) const
{
   const Int_t n = GetNedges() + 1;
   const Int_t nz = GetNz();
   const Int_t nbPnts = nz * n + 2;

   if ((nz < 2) || (nbPnts <= 0) || (n < 2))
      return;

   Int_t c = GetBasicColor();

   Int_t indx = 0, indx1 = 0, indx2 = 0, i, j;

   //  outside circles, number of segments: nz*n
   for (i = 0; i < nz; i++) {
      indx2 = i * n;
      for (j = 1; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j - 1;
         buff.fSegs[indx++] = indx2 + j % (n - 1);
      }
   }

   indx2 = 0;
   // bottom lines
   for (j = 0; j < n; j++) {
      buff.fSegs[indx++] = c;
      buff.fSegs[indx++] = indx2 + j % (n - 1);
      buff.fSegs[indx++] = nbPnts - 2;
   }

   indx2 = (nz - 1) * n;
   // top lines
   for (j = 0; j < n; j++) {
      buff.fSegs[indx++] = c;
      buff.fSegs[indx++] = indx2 + j % (n - 1);
      buff.fSegs[indx++] = nbPnts - 1;
   }

   // outside cylinders, number of segments: (nz-1)*n
   for (i = 0; i < (nz - 1); i++) {
      // outside cylinder
      indx2 = i * n;
      for (j = 0; j < n; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j % (n - 1);
         buff.fSegs[indx++] = indx2 + n + j % (n - 1);
      }
   }

   indx = 0;

   // bottom cap
   indx1 = 0; // start of first z layer
   indx2 = nz * (n - 1);
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 3;
      buff.fPols[indx++] = indx1 + j;
      buff.fPols[indx++] = indx2 + (j + 1) % (n - 1);
      buff.fPols[indx++] = indx2 + j;
   }

   // top cap
   indx1 = (nz - 1) * (n - 1); // start last z layer
   indx2 = nz * (n - 1) + n;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 3;
      buff.fPols[indx++] = indx1 + j; // last z layer
      buff.fPols[indx++] = indx2 + j;
      buff.fPols[indx++] = indx2 + (j + 1) % (n - 1);
   }

   // outside, number of polygons: (nz-1)*(n-1)
   for (Int_t k = 0; k < (nz - 1); k++) {
      indx1 = k * (n - 1);
      indx2 = nz * (n - 1) + n * 2 + k * n;
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = indx1 + j;
         buff.fPols[indx++] = indx2 + j;
         buff.fPols[indx++] = indx1 + j + (n - 1);
         buff.fPols[indx++] = indx2 + (j + 1) % (n - 1);
      }
   }
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Computes projected pgon radius (inner or outer) corresponding to a given Z
// /// value. Fills corresponding coefficients of:
// ///  `Rpg(z) = a + b*z`
// ///
// /// Note: ipl must be in range [0,fNz-2]

// Double_t TGeoVGPgon::Rpg(Double_t z, Int_t ipl, Bool_t inner, Double_t &a, Double_t &b) const
// {
//    Double_t rpg;
//    if (ipl < 0 || ipl > fNz - 2) {
//       Fatal("Rpg", "Plane index parameter ipl=%i out of range\n", ipl);
//       return 0;
//    }
//    Double_t dz = fZ[ipl + 1] - fZ[ipl];
//    if (dz < TGeoShape::Tolerance()) {
//       // radius-changing region
//       rpg = (inner) ? TMath::Min(fRmin[ipl], fRmin[ipl + 1]) : TMath::Max(fRmax[ipl], fRmax[ipl + 1]);
//       a = rpg;
//       b = 0.;
//       return rpg;
//    }
//    Double_t r1 = 0, r2 = 0;
//    if (inner) {
//       r1 = fRmin[ipl];
//       r2 = fRmin[ipl + 1];
//    } else {
//       r1 = fRmax[ipl];
//       r2 = fRmax[ipl + 1];
//    }
//    Double_t dzinv = 1. / dz;
//    a = (r1 * fZ[ipl + 1] - r2 * fZ[ipl]) * dzinv;
//    b = (r2 - r1) * dzinv;
//    return (a + b * z);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Computes projected distance at a given Z for a given ray inside a given sector
// /// and fills coefficients:
// ///   `Rproj = a + b*z`

// Double_t TGeoVGPgon::Rproj(Double_t z, const Double_t *point, const Double_t *dir, Double_t cphi, Double_t sphi,
//                          Double_t &a, Double_t &b) const
// {
//    if (TMath::Abs(dir[2]) < TGeoShape::Tolerance()) {
//       a = b = TGeoShape::Big();
//       return TGeoShape::Big();
//    }
//    Double_t invdirz = 1. / dir[2];
//    a = ((point[0] * dir[2] - point[2] * dir[0]) * cphi + (point[1] * dir[2] - point[2] * dir[1]) * sphi) * invdirz;
//    b = (dir[0] * cphi + dir[1] * sphi) * invdirz;
//    return (a + b * z);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Compute safety from POINT to segment between planes ipl, ipl+1 within safmin.

// Double_t TGeoVGPgon::SafetyToSegment(const Double_t *point, Int_t ipl, Int_t iphi, Bool_t in, Double_t safphi,
//                                    Double_t safmin) const
// {
//    Double_t saf[3];
//    Double_t safe;
//    Int_t i;
//    Double_t r, rpgon, ta, calf;
//    if (ipl < 0 || ipl > fNz - 2)
//       return (safmin + 1.); // error in input plane
//                             // Get info about segment.
//    Double_t dz = fZ[ipl + 1] - fZ[ipl];
//    if (dz < 1E-9)
//       return 1E9; // skip radius-changing segment
//    Double_t znew = point[2] - 0.5 * (fZ[ipl] + fZ[ipl + 1]);
//    saf[0] = 0.5 * dz - TMath::Abs(znew);
//    if (-saf[0] > safmin)
//       return TGeoShape::Big(); // means: stop checking further segments
//    Double_t rmin1 = fRmin[ipl];
//    Double_t rmax1 = fRmax[ipl];
//    Double_t rmin2 = fRmin[ipl + 1];
//    Double_t rmax2 = fRmax[ipl + 1];
//    Double_t divphi = fDphi / fNedges;
//    if (iphi < 0) {
//       Double_t f = 1. / TMath::Cos(0.5 * divphi * TMath::DegToRad());
//       rmax1 *= f;
//       rmax2 *= f;
//       r = TMath::Sqrt(point[0] * point[0] + point[1] * point[1]);
//       Double_t ro1 = 0.5 * (rmin1 + rmin2);
//       Double_t tg1 = (rmin2 - rmin1) / dz;
//       Double_t cr1 = 1. / TMath::Sqrt(1. + tg1 * tg1);
//       Double_t ro2 = 0.5 * (rmax1 + rmax2);
//       Double_t tg2 = (rmax2 - rmax1) / dz;
//       Double_t cr2 = 1. / TMath::Sqrt(1. + tg2 * tg2);
//       Double_t rin = tg1 * znew + ro1;
//       Double_t rout = tg2 * znew + ro2;
//       saf[1] = (ro1 > 0) ? ((r - rin) * cr1) : TGeoShape::Big();
//       saf[2] = (rout - r) * cr2;
//       for (i = 0; i < 3; i++)
//          saf[i] = -saf[i];
//       safe = saf[TMath::LocMax(3, saf)];
//       safe = TMath::Max(safe, safphi);
//       if (safe < 0)
//          safe = 0;
//       return safe;
//    }
//    Double_t ph0 = (fPhi1 + divphi * (iphi + 0.5)) * TMath::DegToRad();
//    r = point[0] * TMath::Cos(ph0) + point[1] * TMath::Sin(ph0);
//    if (rmin1 + rmin2 > 1E-10) {
//       ta = (rmin2 - rmin1) / dz;
//       calf = 1. / TMath::Sqrt(1 + ta * ta);
//       rpgon = rmin1 + (point[2] - fZ[ipl]) * ta;
//       saf[1] = (r - rpgon) * calf;
//    } else {
//       saf[1] = TGeoShape::Big();
//    }
//    ta = (rmax2 - rmax1) / dz;
//    calf = 1. / TMath::Sqrt(1 + ta * ta);
//    rpgon = rmax1 + (point[2] - fZ[ipl]) * ta;
//    saf[2] = (rpgon - r) * calf;
//    if (in) {
//       safe = saf[TMath::LocMin(3, saf)];
//       safe = TMath::Min(safe, safphi);
//    } else {
//       for (i = 0; i < 3; i++)
//          saf[i] = -saf[i];
//       safe = saf[TMath::LocMax(3, saf)];
//       safe = TMath::Max(safe, safphi);
//    }
//    if (safe < 0)
//       safe = 0;
//    return safe;
// }

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGPgon::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   phi1    = " << fPhi1 << ";" << std::endl;
   out << "   dphi    = " << fDphi << ";" << std::endl;
   out << "   nedges = " << fNedges << ";" << std::endl;
   out << "   nz      = " << fNz << ";" << std::endl;
   out << "   auto " << GetPointerName() << " = new TGeoVGPgon(\"" << GetName() << "\", phi1, dphi, nedges, nz);"
       << std::endl;
   for (Int_t i = 0; i < fNz; i++) {
      out << "      z     = " << fZ[i] << ";" << std::endl;
      out << "      rmin  = " << fRmin[i] << ";" << std::endl;
      out << "      rmax  = " << fRmax[i] << ";" << std::endl;
      out << "   " << GetPointerName() << "->DefineSection(" << i << ", z, rmin, rmax);" << std::endl;
   }
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set PGON dimensions starting from an array.

void TGeoVGPgon::SetDimensions(Double_t *param)
{
   fPhi1 = param[0];
   fDphi = param[1];
   fNedges = (Int_t)param[2];
   fNz = (Int_t)param[3];
   if (fNz < 2) {
      Error("SetDimensions", "Pgon %s: Number of Z sections must be > 2", GetName());
      return;
   }
   if (fRmin)
      delete[] fRmin;
   if (fRmax)
      delete[] fRmax;
   if (fZ)
      delete[] fZ;
   fRmin = new Double_t[fNz];
   fRmax = new Double_t[fNz];
   fZ = new Double_t[fNz];
   memset(fRmin, 0, fNz * sizeof(Double_t));
   memset(fRmax, 0, fNz * sizeof(Double_t));
   memset(fZ, 0, fNz * sizeof(Double_t));
   for (Int_t i = 0; i < fNz; i++)
      DefineSection(i, param[4 + 3 * i], param[5 + 3 * i], param[6 + 3 * i]);
}

////////////////////////////////////////////////////////////////////////////////
/// create polygone mesh points

void TGeoVGPgon::SetPoints(Double_t *points) const
{
   Double_t phi, dphi;
   Int_t n = fNedges + 1;
   dphi = fDphi / (n - 1);
   Double_t factor = 1. / TMath::Cos(TMath::DegToRad() * dphi / 2);
   Int_t i, j;
   Int_t indx = 0;

   Bool_t hasInside = HasInsideSurface();

   if (points) {
      for (i = 0; i < GetNz(); i++) {
         if (hasInside)
            for (j = 0; j < n; j++) {
               phi = (fPhi1 + j * dphi) * TMath::DegToRad();
               points[indx++] = factor * fRmin[i] * TMath::Cos(phi);
               points[indx++] = factor * fRmin[i] * TMath::Sin(phi);
               points[indx++] = fZ[i];
            }
         for (j = 0; j < n; j++) {
            phi = (fPhi1 + j * dphi) * TMath::DegToRad();
            points[indx++] = factor * fRmax[i] * TMath::Cos(phi);
            points[indx++] = factor * fRmax[i] * TMath::Sin(phi);
            points[indx++] = fZ[i];
         }
      }

      if (!hasInside) {
         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[0];

         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[GetNz() - 1];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create polygone mesh points

void TGeoVGPgon::SetPoints(Float_t *points) const
{
   Double_t phi, dphi;
   Int_t n = fNedges + 1;
   dphi = fDphi / (n - 1);
   Double_t factor = 1. / TMath::Cos(TMath::DegToRad() * dphi / 2);
   Int_t i, j;
   Int_t indx = 0;

   Bool_t hasInside = HasInsideSurface();

   if (points) {
      for (i = 0; i < fNz; i++) {
         if (hasInside)
            for (j = 0; j < n; j++) {
               phi = (fPhi1 + j * dphi) * TMath::DegToRad();
               points[indx++] = factor * fRmin[i] * TMath::Cos(phi);
               points[indx++] = factor * fRmin[i] * TMath::Sin(phi);
               points[indx++] = fZ[i];
            }
         for (j = 0; j < n; j++) {
            phi = (fPhi1 + j * dphi) * TMath::DegToRad();
            points[indx++] = factor * fRmax[i] * TMath::Cos(phi);
            points[indx++] = factor * fRmax[i] * TMath::Sin(phi);
            points[indx++] = fZ[i];
         }
      }

      if (!hasInside) {
         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[0];

         points[indx++] = 0;
         points[indx++] = 0;
         points[indx++] = fZ[GetNz() - 1];
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGPgon::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   nvert = nsegs = npols = 0;

   Int_t n = GetNedges() + 1;
   Int_t nz = GetNz();

   if (nz < 2)
      return;

   if (HasInsideSurface()) {
      Bool_t specialCase = TGeoShape::IsSameWithinTolerance(GetDphi(), 360);
      nvert = nz * 2 * n;
      nsegs = 4 * (nz * n - 1 + (specialCase ? 1 : 0));
      npols = 2 * (nz * n - 1 + (specialCase ? 1 : 0));
   } else {
      nvert = nz * n + 2;
      nsegs = nz * (n - 1) + n * 2 + (nz - 1) * n;
      npols = 2 * (n - 1) + (nz - 1) * (n - 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGPgon::GetNmeshVertices() const
{
   Int_t nvert, nsegs, npols;

   GetMeshNumbers(nvert, nsegs, npols);

   return nvert;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGPgon::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Returns true when pgon has internal surface
/// It will be only disabled when all Rmin values are 0 
/// From TGeoPcon

Bool_t TGeoVGPgon::HasInsideSurface() const
{
   // only when full 360 is used, internal part can be excluded
   Bool_t specialCase = TGeoShape::IsSameWithinTolerance(GetDphi(), 360);
   if (!specialCase)
      return kTRUE;

   for (Int_t i = 0; i < GetNz(); i++)
      if (fRmin[i] > 0.)
         return kTRUE;

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGPgon::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t nbPnts, nbSegs, nbPols;
      GetMeshNumbers(nbPnts, nbSegs, nbPols);
      if (nbPnts > 0) {
         if (buffer.SetRawSizes(nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols)) {
            buffer.SetSectionsValid(TBuffer3D::kRawSizes);
         }
      }
   }
   // TODO: Push down to TGeoShape?? Would have to do raw sizes set first..
   // can rest of TGeoShape be deferred until after this?
   if ((reqSections & TBuffer3D::kRaw) && buffer.SectionsValid(TBuffer3D::kRawSizes)) {
      SetPoints(buffer.fPnts);
      if (!buffer.fLocalFrame) {
         TransformPoints(buffer.fPnts, buffer.NbPnts());
      }

      SetSegsAndPols(buffer);
      buffer.SetSectionsValid(TBuffer3D::kRaw);
   }

   return buffer;
}

#endif // ROOT_USE_VECGEOM_SOLIDS