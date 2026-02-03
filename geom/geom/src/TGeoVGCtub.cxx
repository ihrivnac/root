// @(#)root/geom:$Id$
// Author: Andrei Gheata   24/10/01
// TGeoTube::Contains() and DistFromInside/In() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoVGCtub.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGCtub);

// TGeoVGCtub::TGeoVGCtub()
// {
//    // default ctor
//    fNlow[0] = fNlow[1] = fNhigh[0] = fNhigh[1] = 0.;
//    fNlow[2] = -1;
//    fNhigh[2] = 1;
// }

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoVGCtub::TGeoVGCtub(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2, Double_t lx, Double_t ly,
                   Double_t lz, Double_t tx, Double_t ty, Double_t tz)
   : Base_t("", rmin, rmax, dz, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad(), lx, ly, lz, tx, ty, tz)
{
   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   SetShapeBit(kGeoCtub);
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor

TGeoVGCtub::TGeoVGCtub(const char *name, Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2,
                   Double_t lx, Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
   : Base_t(name, rmin, rmax, dz, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad(), lx, ly, lz, tx, ty, tz)
{
   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
   SetShapeBit(kGeoCtub);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// ctor with parameters

TGeoVGCtub::TGeoVGCtub(Double_t *params) : TGeoTubeSeg(0, 0, 0, 0, 0)
{
   SetCtubDimensions(params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8],
                     params[9], params[10]);
   SetShapeBit(kGeoCtub);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGCtub::~TGeoVGCtub() {}

////////////////////////////////////////////////////////////////////////////////
/// Computes capacity of the shape in [length^3]

Double_t TGeoVGCtub::Capacity() const
{
   Double_t capacity = TGeoTubeSeg::Capacity();
   return capacity;
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Init frequently used trigonometric values

void TGeoVGTubeSeg::InitTrigonometry()
{
   Double_t phi1 = sphi();
   Double_t phi2 = (sphi() + dphi());
   fC1 = TMath::Cos(phi1);
   fS1 = TMath::Sin(phi1);
   fC2 = TMath::Cos(phi2);
   fS2 = TMath::Sin(phi2);
   Double_t fio = 0.5 * (phi1 + phi2);
   fCm = TMath::Cos(fio);
   fSm = TMath::Sin(fio);
   Double_t dfi = 0.5 * (phi2 - phi1);
   fCdfi = TMath::Cos(dfi);
}

////////////////////////////////////////////////////////////////////////////////
/// compute minimum bounding box of the ctub

void TGeoVGCtub::ComputeBBox()
{
   TGeoTubeSeg::ComputeBBox();
   if ((fNlow[2] > -(1E-10)) || (fNhigh[2] < 1E-10)) {
      Error("ComputeBBox", "In shape %s wrong definition of cut planes", GetName());
      return;
   }
   Double_t xc = 0, yc = 0;
   Double_t zmin = 0, zmax = 0;
   Double_t z1;
   Double_t z[8];
   // check if nxy is in the phi range
   Double_t phi_low = TMath::ATan2(fNlow[1], fNlow[0]) * TMath::RadToDeg();
   Double_t phi_hi = TMath::ATan2(fNhigh[1], fNhigh[0]) * TMath::RadToDeg();
   Bool_t in_range_low = kFALSE;
   Bool_t in_range_hi = kFALSE;

   Int_t i;
   for (i = 0; i < 2; i++) {
      if (phi_low < 0)
         phi_low += 360.;
      Double_t dphi = fPhi2 - fPhi1;
      if (dphi < 0)
         dphi += 360.;
      Double_t ddp = phi_low - fPhi1;
      if (ddp < 0)
         ddp += 360.;
      if (ddp <= dphi) {
         xc = rmin() * TMath::Cos(phi_low * TMath::DegToRad());
         yc = rmin() * TMath::Sin(phi_low * TMath::DegToRad());
         z1 = GetZcoord(xc, yc, -z());
         xc = rmax() * TMath::Cos(phi_low * TMath::DegToRad());
         yc = rmax() * TMath::Sin(phi_low * TMath::DegToRad());
         z1 = TMath::Min(z1, GetZcoord(xc, yc, -z()));
         if (in_range_low)
            zmin = TMath::Min(zmin, z1);
         else
            zmin = z1;
         in_range_low = kTRUE;
      }
      phi_low += 180;
      if (phi_low > 360)
         phi_low -= 360.;
   }

   for (i = 0; i < 2; i++) {
      if (phi_hi < 0)
         phi_hi += 360.;
      Double_t dphi = fPhi2 - fPhi1;
      if (dphi < 0)
         dphi += 360.;
      Double_t ddp = phi_hi - fPhi1;
      if (ddp < 0)
         ddp += 360.;
      if (ddp <= dphi) {
         xc = rmin() * TMath::Cos(phi_hi * TMath::DegToRad());
         yc = rmin() * TMath::Sin(phi_hi * TMath::DegToRad());
         z1 = GetZcoord(xc, yc, z());
         xc = rmax() * TMath::Cos(phi_hi * TMath::DegToRad());
         yc = rmax() * TMath::Sin(phi_hi * TMath::DegToRad());
         z1 = TMath::Max(z1, GetZcoord(xc, yc, z()));
         if (in_range_hi)
            zmax = TMath::Max(zmax, z1);
         else
            zmax = z1;
         in_range_hi = kTRUE;
      }
      phi_hi += 180;
      if (phi_hi > 360)
         phi_hi -= 360.;
   }

   xc = rmin() * fC1;
   yc = rmin() * fS1;
   z[0] = GetZcoord(xc, yc, -z());
   z[4] = GetZcoord(xc, yc, z());

   xc = rmin() * fC2;
   yc = rmin() * fS2;
   z[1] = GetZcoord(xc, yc, -z());
   z[5] = GetZcoord(xc, yc, z());

   xc = rmax() * fC1;
   yc = rmax() * fS1;
   z[2] = GetZcoord(xc, yc, -z());
   z[6] = GetZcoord(xc, yc, z());

   xc = rmax() * fC2;
   yc = rmax() * fS2;
   z[3] = GetZcoord(xc, yc, -z());
   z[7] = GetZcoord(xc, yc, z());

   z1 = z[TMath::LocMin(4, &z[0])];
   if (in_range_low)
      zmin = TMath::Min(zmin, z1);
   else
      zmin = z1;

   z1 = z[TMath::LocMax(4, &z[4]) + 4];
   if (in_range_hi)
      zmax = TMath::Max(zmax, z1);
   else
      zmax = z1;

   fDZ = 0.5 * (zmax - zmin);
   fOrigin[2] = 0.5 * (zmax + zmin);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute normal to closest surface from POINT.

// void TGeoVGCtub::ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const
// {
//    Double_t saf[4];
//    Bool_t isseg = kTRUE;
//    if (TMath::Abs(fPhi2 - fPhi1 - 360.) < 1E-8)
//       isseg = kFALSE;
//    Double_t rsq = point[0] * point[0] + point[1] * point[1];
//    Double_t r = TMath::Sqrt(rsq);

//    saf[0] = TMath::Abs(point[0] * fNlow[0] + point[1] * fNlow[1] + (z() + point[2]) * fNlow[2]);
//    saf[1] = TMath::Abs(point[0] * fNhigh[0] + point[1] * fNhigh[1] - (z() - point[2]) * fNhigh[2]);
//    saf[2] = (rmin() > 1E-10) ? TMath::Abs(r - rmin()) : TGeoShape::Big();
//    saf[3] = TMath::Abs(rmax() - r);
//    Int_t i = TMath::LocMin(4, saf);
//    if (isseg) {
//       if (TGeoShape::IsCloseToPhi(saf[i], point, fC1, fS1, fC2, fS2)) {
//          TGeoShape::NormalPhi(point, dir, norm, fC1, fS1, fC2, fS2);
//          return;
//       }
//    }
//    if (i == 0) {
//       memcpy(norm, fNlow, 3 * sizeof(Double_t));
//       if (norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2] < 0) {
//          norm[0] = -norm[0];
//          norm[1] = -norm[1];
//          norm[2] = -norm[2];
//       }
//       return;
//    }
//    if (i == 1) {
//       memcpy(norm, fNhigh, 3 * sizeof(Double_t));
//       if (norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2] < 0) {
//          norm[0] = -norm[0];
//          norm[1] = -norm[1];
//          norm[2] = -norm[2];
//       }
//       return;
//    }

//    norm[2] = 0;
//    Double_t phi = TMath::ATan2(point[1], point[0]);
//    norm[0] = TMath::Cos(phi);
//    norm[1] = TMath::Sin(phi);
//    if (norm[0] * dir[0] + norm[1] * dir[1] < 0) {
//       norm[0] = -norm[0];
//       norm[1] = -norm[1];
//    }
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// check if point is contained in the cut tube
// /// check the lower cut plane

// Bool_t TGeoVGCtub::Contains(const Double_t *point) const
// {
//    Double_t zin = point[0] * fNlow[0] + point[1] * fNlow[1] + (point[2] + z()) * fNlow[2];
//    if (zin > 0)
//       return kFALSE;
//    // check the higher cut plane
//    zin = point[0] * fNhigh[0] + point[1] * fNhigh[1] + (point[2] - z()) * fNhigh[2];
//    if (zin > 0)
//       return kFALSE;
//    // check radius
//    Double_t r2 = point[0] * point[0] + point[1] * point[1];
//    if ((r2 < rmin() * rmin()) || (r2 > rmax() * rmax()))
//       return kFALSE;
//    // check phi
//    Double_t phi = TMath::ATan2(point[1], point[0]) * TMath::RadToDeg();
//    if (phi < 0)
//       phi += 360.;
//    Double_t dphi = fPhi2 - fPhi1;
//    Double_t ddp = phi - fPhi1;
//    if (ddp < 0)
//       ddp += 360.;
//    //   if (ddp>360) ddp-=360;
//    if (ddp > dphi)
//       return kFALSE;
//    return kTRUE;
// }

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGCtub::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 1:
      xlo = rmin();
      xhi = rmax();
      dx = xhi - xlo;
      return dx;
   case 2:
      xlo = fPhi1;
      xhi = fPhi2;
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// compute real Z coordinate of a point belonging to either lower or
/// higher caps (z should be either +z() or -z())

Double_t TGeoVGCtub::GetZcoord(Double_t xc, Double_t yc, Double_t zc) const
{
   Double_t newz = 0;
   if (zc < 0)
      newz = -z() - (xc * fNlow[0] + yc * fNlow[1]) / fNlow[2];
   else
      newz = z() - (xc * fNhigh[0] + yc * fNhigh[1]) / fNhigh[2];
   return newz;
}

////////////////////////////////////////////////////////////////////////////////
/// compute distance from outside point to surface of the cut tube

// Double_t
// TGeoVGCtub::DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
// {
//    if (iact < 3 && safe) {
//       *safe = Safety(point, kFALSE);
//       if (iact == 0)
//          return TGeoShape::Big();
//       if ((iact == 1) && (step <= *safe))
//          return TGeoShape::Big();
//    }
//    // Check if the bounding box is crossed within the requested distance
//    Double_t sdist = TGeoBBox::DistFromOutside(point, dir, fDX, fDY, fDZ, fOrigin, step);
//    if (sdist >= step)
//       return TGeoShape::Big();
//    Double_t saf[2];
//    saf[0] = point[0] * fNlow[0] + point[1] * fNlow[1] + (z() + point[2]) * fNlow[2];
//    saf[1] = point[0] * fNhigh[0] + point[1] * fNhigh[1] + (point[2] - z()) * fNhigh[2];
//    Double_t rsq = point[0] * point[0] + point[1] * point[1];
//    Double_t r = TMath::Sqrt(rsq);
//    Double_t cpsi = 0;
//    Bool_t tub = kFALSE;
//    if (TMath::Abs(fPhi2 - fPhi1 - 360.) < 1E-8)
//       tub = kTRUE;

//    // find distance to shape
//    Double_t r2;
//    Double_t calf = dir[0] * fNlow[0] + dir[1] * fNlow[1] + dir[2] * fNlow[2];
//    // check Z planes
//    Double_t xi, yi, zi, s;
//    if (saf[0] > 0) {
//       if (calf < 0) {
//          s = -saf[0] / calf;
//          xi = point[0] + s * dir[0];
//          yi = point[1] + s * dir[1];
//          r2 = xi * xi + yi * yi;
//          if (((rmin() * rmin()) <= r2) && (r2 <= (rmax() * rmax()))) {
//             if (tub)
//                return s;
//             cpsi = (xi * fCm + yi * fSm) / TMath::Sqrt(r2);
//             if (cpsi >= fCdfi)
//                return s;
//          }
//       }
//    }
//    calf = dir[0] * fNhigh[0] + dir[1] * fNhigh[1] + dir[2] * fNhigh[2];
//    if (saf[1] > 0) {
//       if (calf < 0) {
//          s = -saf[1] / calf;
//          xi = point[0] + s * dir[0];
//          yi = point[1] + s * dir[1];
//          r2 = xi * xi + yi * yi;
//          if (((rmin() * rmin()) <= r2) && (r2 <= (rmax() * rmax()))) {
//             if (tub)
//                return s;
//             cpsi = (xi * fCm + yi * fSm) / TMath::Sqrt(r2);
//             if (cpsi >= fCdfi)
//                return s;
//          }
//       }
//    }

//    // check outer cyl. surface
//    Double_t nsq = dir[0] * dir[0] + dir[1] * dir[1];
//    if (TMath::Abs(nsq) < 1E-10)
//       return TGeoShape::Big();
//    Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
//    Double_t b, d;
//    // only r>rmax() coming inwards has to be considered
//    if (r > rmax() && rdotn < 0) {
//       TGeoTube::DistToTube(rsq, nsq, rdotn, rmax(), b, d);
//       if (d > 0) {
//          s = -b - d;
//          if (s > 0) {
//             xi = point[0] + s * dir[0];
//             yi = point[1] + s * dir[1];
//             zi = point[2] + s * dir[2];
//             if ((-xi * fNlow[0] - yi * fNlow[1] - (zi + z()) * fNlow[2]) > 0) {
//                if ((-xi * fNhigh[0] - yi * fNhigh[1] + (z() - zi) * fNhigh[2]) > 0) {
//                   if (tub)
//                      return s;
//                   cpsi = (xi * fCm + yi * fSm) / rmax();
//                   if (cpsi >= fCdfi)
//                      return s;
//                }
//             }
//          }
//       }
//    }
//    // check inner cylinder
//    Double_t snxt = TGeoShape::Big();
//    if (rmin() > 0) {
//       TGeoTube::DistToTube(rsq, nsq, rdotn, rmin(), b, d);
//       if (d > 0) {
//          s = -b + d;
//          if (s > 0) {
//             xi = point[0] + s * dir[0];
//             yi = point[1] + s * dir[1];
//             zi = point[2] + s * dir[2];
//             if ((-xi * fNlow[0] - yi * fNlow[1] - (zi + z()) * fNlow[2]) > 0) {
//                if ((-xi * fNhigh[0] - yi * fNhigh[1] + (z() - zi) * fNhigh[2]) > 0) {
//                   if (tub)
//                      return s;
//                   cpsi = (xi * fCm + yi * fSm) / rmin();
//                   if (cpsi >= fCdfi)
//                      snxt = s;
//                }
//             }
//          }
//       }
//    }
//    // check phi planes
//    if (tub)
//       return snxt;
//    Double_t un = dir[0] * fS1 - dir[1] * fC1;
//    if (un < -TGeoShape::Tolerance()) {
//       s = (point[1] * fC1 - point[0] * fS1) / un;
//       if (s >= 0) {
//          xi = point[0] + s * dir[0];
//          yi = point[1] + s * dir[1];
//          zi = point[2] + s * dir[2];
//          if ((-xi * fNlow[0] - yi * fNlow[1] - (zi + z()) * fNlow[2]) > 0) {
//             if ((-xi * fNhigh[0] - yi * fNhigh[1] + (z() - zi) * fNhigh[2]) > 0) {
//                r2 = xi * xi + yi * yi;
//                if ((rmin() * rmin() <= r2) && (r2 <= rmax() * rmax())) {
//                   if ((yi * fCm - xi * fSm) <= 0) {
//                      if (s < snxt)
//                         snxt = s;
//                   }
//                }
//             }
//          }
//       }
//    }
//    un = dir[0] * fS2 - dir[1] * fC2;
//    if (un > TGeoShape::Tolerance()) {
//       s = (point[1] * fC2 - point[0] * fS2) / un;
//       if (s >= 0) {
//          xi = point[0] + s * dir[0];
//          yi = point[1] + s * dir[1];
//          zi = point[2] + s * dir[2];
//          if ((-xi * fNlow[0] - yi * fNlow[1] - (zi + z()) * fNlow[2]) > 0) {
//             if ((-xi * fNhigh[0] - yi * fNhigh[1] + (z() - zi) * fNhigh[2]) > 0) {
//                r2 = xi * xi + yi * yi;
//                if ((rmin() * rmin() <= r2) && (r2 <= rmax() * rmax())) {
//                   if ((yi * fCm - xi * fSm) >= 0) {
//                      if (s < snxt)
//                         snxt = s;
//                   }
//                }
//             }
//          }
//       }
//    }
//    return snxt;
// }

////////////////////////////////////////////////////////////////////////////////
/// compute distance from inside point to surface of the cut tube

// Double_t
// TGeoVGCtub::DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact, Double_t step, Double_t *safe) const
// {
//    if (iact < 3 && safe)
//       *safe = Safety(point, kTRUE);
//    if (iact == 0)
//       return TGeoShape::Big();
//    if ((iact == 1) && (*safe > step))
//       return TGeoShape::Big();
//    Double_t rsq = point[0] * point[0] + point[1] * point[1];
//    Bool_t tub = kFALSE;
//    if (TMath::Abs(fPhi2 - fPhi1 - 360.) < 1E-8)
//       tub = kTRUE;
//    // compute distance to surface
//    // Do Z
//    Double_t sz = TGeoShape::Big();
//    Double_t saf[2];
//    saf[0] = -point[0] * fNlow[0] - point[1] * fNlow[1] - (z() + point[2]) * fNlow[2];
//    saf[1] = -point[0] * fNhigh[0] - point[1] * fNhigh[1] + (z() - point[2]) * fNhigh[2];
//    Double_t calf = dir[0] * fNlow[0] + dir[1] * fNlow[1] + dir[2] * fNlow[2];
//    if (calf > 0)
//       sz = saf[0] / calf;

//    calf = dir[0] * fNhigh[0] + dir[1] * fNhigh[1] + dir[2] * fNhigh[2];
//    if (calf > 0) {
//       Double_t sz1 = saf[1] / calf;
//       if (sz1 < sz)
//          sz = sz1;
//    }

//    // Do R
//    Double_t nsq = dir[0] * dir[0] + dir[1] * dir[1];
//    // track parallel to Z
//    if (TMath::Abs(nsq) < 1E-10)
//       return sz;
//    Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
//    Double_t sr = TGeoShape::Big();
//    Double_t b, d;
//    Bool_t skip_outer = kFALSE;
//    // inner cylinder
//    if (rmin() > 1E-10) {
//       TGeoTube::DistToTube(rsq, nsq, rdotn, rmin(), b, d);
//       if (d > 0) {
//          sr = -b - d;
//          if (sr > 0)
//             skip_outer = kTRUE;
//       }
//    }
//    // outer cylinder
//    if (!skip_outer) {
//       TGeoTube::DistToTube(rsq, nsq, rdotn, rmax(), b, d);
//       if (d > 0) {
//          sr = -b + d;
//          if (sr < 0)
//             sr = TGeoShape::Big();
//       } else {
//          return 0.; // already outside
//       }
//    }
//    // phi planes
//    Double_t sfmin = TGeoShape::Big();
//    if (!tub)
//       sfmin = TGeoShape::DistToPhiMin(point, dir, fS1, fC1, fS2, fC2, fSm, fCm);
//    return TMath::Min(TMath::Min(sz, sr), sfmin);
// }

////////////////////////////////////////////////////////////////////////////////
/// Divide the tube along one axis.

TGeoVolume *TGeoVGCtub::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Warning("Divide", "In shape %s division of a cut tube not implemented", GetName());
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGCtub::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (!mother->TestShapeBit(kGeoTube)) {
      Error("GetMakeRuntimeShape", "Invalid mother for shape %s", GetName());
      return nullptr;
   }
   Double_t rmin, rmax, dz;
   rmin = Base_t::rmin();
   rmax = Base_t::rmax();
   dz = z();
   if (z() < 0)
      dz = ((TGeoTVGCtub *)mother)->GetDz();
   if (rmin() < 0)
      rmin = ((TGeoVGCtub *)mother)->GetRmin();
   if ((rmax() < 0) || (rmax() <= rmin()))
      rmax = ((TGeoVGCtub *)mother)->GetRmax();

   return (new TGeoVGCtub(rmin, rmax, dz, fPhi1, fPhi2, fNlow[0], fNlow[1], fNlow[2], fNhigh[0], fNhigh[1], fNhigh[2]));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGCtub::InspectShape() const
{
   printf("*** Shape %s: TGeoVGCtub ***\n", GetName());
   printf("    Rmin = %11.5f\n", rmin());
   printf("    Rmax = %11.5f\n", rmax());
   printf("    dz   = %11.5f\n", z());
   printf("    phi1 = %11.5f\n", sphi()* TMath::RadToDeg());
   printf("    phi2 = %11.5f\n", (sphi() + dphi()) * TMath::RadToDeg());
   printf("    lx = %11.5f\n", fNlow[0]);
   printf("    ly = %11.5f\n", fNlow[1]);
   printf("    lz = %11.5f\n", fNlow[2]);
   printf("    tx = %11.5f\n", fNhigh[0]);
   printf("    ty = %11.5f\n", fNhigh[1]);
   printf("    tz = %11.5f\n", fNhigh[2]);
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

// ////////////////////////////////////////////////////////////////////////////////
// /// computes the closest distance from given point to this shape, according
// /// to option. The matching point on the shape is stored in spoint.

// Double_t TGeoVGCtub::Safety(const Double_t *point, Bool_t in) const
// {
//    Double_t saf[4];
//    Double_t rsq = point[0] * point[0] + point[1] * point[1];
//    Double_t r = TMath::Sqrt(rsq);
//    Bool_t isseg = kTRUE;
//    if (TMath::Abs(fPhi2 - fPhi1 - 360.) < 1E-8)
//       isseg = kFALSE;

//    saf[0] = -point[0] * fNlow[0] - point[1] * fNlow[1] - (z() + point[2]) * fNlow[2];
//    saf[1] = -point[0] * fNhigh[0] - point[1] * fNhigh[1] + (z() - point[2]) * fNhigh[2];
//    saf[2] = (rmin() < 1E-10 && !isseg) ? TGeoShape::Big() : (r - rmin());
//    saf[3] = rmax() - r;
//    Double_t safphi = TGeoShape::Big();
//    if (isseg)
//       safphi = TGeoShape::SafetyPhi(point, in, fPhi1, fPhi2);

//    if (in) {
//       Double_t safe = saf[TMath::LocMin(4, saf)];
//       return TMath::Min(safe, safphi);
//    }
//    for (Int_t i = 0; i < 4; i++)
//       saf[i] = -saf[i];
//    Double_t safe = saf[TMath::LocMax(4, saf)];
//    if (isseg)
//       return TMath::Max(safe, safphi);
//    return safe;
// }

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGCtub::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << rmin() << ";" << std::endl;
   out << "   rmax = " << rmax() << ";" << std::endl;
   out << "   dz   = " << z() << ";" << std::endl;
   out << "   phi1 = " << sphi() * TMath::RadToDeg() << ";" << std::endl;
   out << "   phi2 = " << (sphi() + dphi()) * TMath::RadToDeg() << ";" << std::endl;
   out << "   lx   = " << fNlow[0] << ";" << std::endl;
   out << "   ly   = " << fNlow[1] << ";" << std::endl;
   out << "   lz   = " << fNlow[2] << ";" << std::endl;
   out << "   tx   = " << fNhigh[0] << ";" << std::endl;
   out << "   ty   = " << fNhigh[1] << ";" << std::endl;
   out << "   tz   = " << fNhigh[2] << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGCtub(\"" << GetName()
       << "\",rmin,rmax,dz,phi1,phi2,lx,ly,lz,tx,ty,tz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi1

Double_t TGeoVGCtub::GetPhi1() const
{
   return sphi() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi2

Double_t TGeoVGCtub::GetPhi2() const
{
   return (sphi() + dphi()) * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the tube segment.
/// The segment will be from phiStart to phiEnd expressed in degree.

void TGeoVGTubeSeg::SetTubsDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart, Double_t phiEnd)
{
   SetRMin(rmin);
   SetRMin(rmax);
   SetDz(dz);
   auto phi1 = phiStart;
   if (phi1 < 0)
      phi1 += 360.;
   SetSPhi(phi1 * TMath::DegToRad());
   auto phi2 = phiEnd;
   while (phi2 <= phi1)
      phi2 += 360.;
   SetDPhi((phi2 - phi1) * TMath::DegToRad());
   if (TGeoShape::IsSameWithinTolerance(sphi()* TMath::RadToDeg(), (sphi() + dphi()) * TMath::RadToDeg()))
      Fatal("SetTubsDimensions", "In shape %s invalid phi1=%g, phi2=%g\n", GetName(), phi1, phi2);
   InitTrigonometry();
}

////////////////////////////////////////////////////////////////////////////////
/// set dimensions of a cut tube

void TGeoVGCtub::SetCtubDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phi1, Double_t phi2, Double_t lx,
                                 Double_t ly, Double_t lz, Double_t tx, Double_t ty, Double_t tz)
{
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);

   fNlow[0] = lx;
   fNlow[1] = ly;
   fNlow[2] = lz;
   fNhigh[0] = tx;
   fNhigh[1] = ty;
   fNhigh[2] = tz;
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the cut tube starting from a list.

void TGeoVGCtub::SetDimensions(Double_t *param)
{
   SetCtubDimensions(param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9],
                     param[10]);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoVGCtub::GetPointsOnSegments(Int_t /*npoints*/, Double_t * /*array*/) const
{
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create mesh points for the cut tube.

void TGeoVGCtub::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2 < phi1)
      phi2 += 360.;
   n = gGeoManager->GetNsegments() + 1;

   dphi = (phi2 - phi1) / (n - 1);
   dz = z();

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi) * TMath::DegToRad();
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = GetZcoord(points[indx - 2], points[indx - 1], dz);
         points[indx] = GetZcoord(points[indx - 2], points[indx - 1], -dz);
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi) * TMath::DegToRad();
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = GetZcoord(points[indx - 2], points[indx - 1], dz);
         points[indx] = GetZcoord(points[indx - 2], points[indx - 1], -dz);
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create mesh points for the cut tube.

void TGeoVGCtub::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = fPhi1;
   phi2 = fPhi2;
   if (phi2 < phi1)
      phi2 += 360.;
   n = gGeoManager->GetNsegments() + 1;

   dphi = (phi2 - phi1) / (n - 1);
   dz = z();

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi) * TMath::DegToRad();
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = GetZcoord(points[indx - 2], points[indx - 1], dz);
         points[indx] = GetZcoord(points[indx - 2], points[indx - 1], -dz);
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi) * TMath::DegToRad();
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = GetZcoord(points[indx - 2], points[indx - 1], dz);
         points[indx] = GetZcoord(points[indx - 2], points[indx - 1], -dz);
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGCtub::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   // TGeoTubeSeg::GetMeshNumbers(nvert, nsegs, npols);
   Int_t n = gGeoManager->GetNsegments() + 1;
   nvert = n * 4;
   nsegs = n * 8;
   npols = n * 4 - 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGCtub::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t numPoints = n * 4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGCtub::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DCutTube buffer;

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      // These from TBuffer3DCutTube / TGeoVGCtub
      buffer.fRadiusInner = rmin();
      buffer.fRadiusOuter = rmax();
      buffer.fHalfLength = z();
      buffer.fPhiMin = fPhi1;
      buffer.fPhiMax = fPhi2;

      for (UInt_t i = 0; i < 3; i++) {
         buffer.fLowPlaneNorm[i] = fNlow[i];
         buffer.fHighPlaneNorm[i] = fNhigh[i];
      }
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments() + 1;
      Int_t nbPnts = 4 * n;
      Int_t nbSegs = 2 * nbPnts;
      Int_t nbPols = nbPnts - 2;
      if (buffer.SetRawSizes(nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
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

// ////////////////////////////////////////////////////////////////////////////////
// /// Check the inside status for each of the points in the array.
// /// Input: Array of point coordinates + vector size
// /// Output: Array of Booleans for the inside of each point

// void TGeoVGCtub::Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const
// {
//    for (Int_t i = 0; i < vecsize; i++)
//       inside[i] = Contains(&points[3 * i]);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Compute the normal for an array o points so that norm.dot.dir is positive
// /// Input: Arrays of point coordinates and directions + vector size
// /// Output: Array of normal directions

// void TGeoVGCtub::ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize)
// {
//    for (Int_t i = 0; i < vecsize; i++)
//       ComputeNormal(&points[3 * i], &dirs[3 * i], &norms[3 * i]);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Compute distance from array of input points having directions specified by dirs. Store output in dists

// void TGeoVGCtub::DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
//                                 Double_t *step) const
// {
//    for (Int_t i = 0; i < vecsize; i++)
//       dists[i] = DistFromInside(&points[3 * i], &dirs[3 * i], 3, step[i]);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Compute distance from array of input points having directions specified by dirs. Store output in dists

// void TGeoVGCtub::DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
//                                  Double_t *step) const
// {
//    for (Int_t i = 0; i < vecsize; i++)
//       dists[i] = DistFromOutside(&points[3 * i], &dirs[3 * i], 3, step[i]);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// Compute safe distance from each of the points in the input array.
// /// Input: Array of point coordinates, array of statuses for these points, size of the arrays
// /// Output: Safety values

// void TGeoVGCtub::Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const
// {
//    for (Int_t i = 0; i < vecsize; i++)
//       safe[i] = Safety(&points[3 * i], inside[i]);
// }
