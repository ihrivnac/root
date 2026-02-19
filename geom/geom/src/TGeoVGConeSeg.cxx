// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoVGConeSeg::Contains() and DistFromInside() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TGeoCone.h"
#include "TGeoVGCone.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

ClassImp(TGeoVGConeSeg);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGConeSeg::TGeoVGConeSeg()
 : Base_t("", 0., 0., 0., 0., 0., 0., vecgeom::kTwoPi),
   fS1(0.),
   fC1(0.),
   fS2(0.),
   fC2(0.),
   fSm(0.),
   fCm(0.),
   fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoConeSeg);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGConeSeg::TGeoVGConeSeg(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t phi1,
                         Double_t phi2)
 : Base_t("", rmin1, rmax1, rmin2, rmax2, dz, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad()),
   fS1(0.),
   fC1(0.),
   fS2(0.),
   fC2(0.),
   fSm(0.),
   fCm(0.),
   fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoConeSeg);
   SetConsDimensions(dz, rmin1, rmax1, rmin2, rmax2, phi1, phi2);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGConeSeg::TGeoVGConeSeg(const char *name, Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2,
                         Double_t phi1, Double_t phi2)
 : Base_t(name, rmin1, rmax1, rmin2, rmax2, dz, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad()),
   fS1(0.),
   fC1(0.),
   fS2(0.),
   fC2(0.),
   fSm(0.),
   fCm(0.),
   fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoConeSeg);
   SetConsDimensions(dz, rmin1, rmax1, rmin2, rmax2, phi1, phi2);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = dz
///  - param[1] = Rmin1
///  - param[2] = Rmax1
///  - param[3] = Rmin2
///  - param[4] = Rmax2
///  - param[5] = phi1
///  - param[6] = phi2

TGeoVGConeSeg::TGeoVGConeSeg(Double_t *param)
 : Base_t("", param[1], param[2], param[3], param[3], param[0], param[5] * TMath::DegToRad(), (param[6] - param[5]) * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoConeSeg);
   SetDimensions(param);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGConeSeg::~TGeoVGConeSeg() {}

////////////////////////////////////////////////////////////////////////////////
/// Init frequently used trigonometric values

void TGeoVGConeSeg::InitTrigonometry()
{
   Double_t phi1 = GetSPhi();
   Double_t phi2 = GetSPhi() + GetDPhi();
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
/// Static: Computes capacity of the shape in [length^3]

Double_t TGeoVGConeSeg::Capacity(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2,
                               Double_t phi1, Double_t phi2)
{
   Double_t capacity = (TMath::Abs(phi2 - phi1) * TMath::DegToRad() * dz / 3.) *
                       (rmax1 * rmax1 + rmax2 * rmax2 + rmax1 * rmax2 - rmin1 * rmin1 - rmin2 * rmin2 - rmin1 * rmin2);
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the tube segment

void TGeoVGConeSeg::ComputeBBox()
{
   Double_t rmin, rmax;
   rmin = TMath::Min(GetRmin1(), GetRmin2());
   rmax = TMath::Max(GetRmax1(), GetRmax2());

   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax * fC1;
   yc[0] = rmax * fS1;
   xc[1] = rmax * fC2;
   yc[1] = rmax * fS2;
   xc[2] = rmin * fC1;
   yc[2] = rmin * fS1;
   xc[3] = rmin * fC2;
   yc[3] = rmin * fS2;

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];

   Double_t dp = GetDPhi() * TMath::RadToDeg();
   Double_t ddp = -GetSPhi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp <= dp)
      xmax = rmax;
   ddp = 90 - GetSPhi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp <= dp)
      ymax = rmax;
   ddp = 180 - GetSPhi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp <= dp)
      xmin = -rmax;
   ddp = 270 - GetSPhi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp <= dp)
      ymin = -rmax;
   fOrigin[0] = (xmax + xmin) / 2;
   fOrigin[1] = (ymax + ymin) / 2;
   fOrigin[2] = 0;
   fDX = (xmax - xmin) / 2;
   fDY = (ymax - ymin) / 2;
   fDZ = GetDz();
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute normal to closest surface from POINT.

void TGeoVGConeSeg::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t dz,
                                 Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1,
                                 Double_t s1, Double_t c2, Double_t s2)
{
   Double_t saf[2];
   Double_t ro1 = 0.5 * (rmin1 + rmin2);
   Double_t tg1 = 0.5 * (rmin2 - rmin1) / dz;
   Double_t cr1 = 1. / TMath::Sqrt(1. + tg1 * tg1);
   Double_t ro2 = 0.5 * (rmax1 + rmax2);
   Double_t tg2 = 0.5 * (rmax2 - rmax1) / dz;
   Double_t cr2 = 1. / TMath::Sqrt(1. + tg2 * tg2);

   Double_t r = TMath::Sqrt(point[0] * point[0] + point[1] * point[1]);
   Double_t rin = tg1 * point[2] + ro1;
   Double_t rout = tg2 * point[2] + ro2;
   saf[0] = (ro1 > 0) ? (TMath::Abs((r - rin) * cr1)) : TGeoShape::Big();
   saf[1] = TMath::Abs((rout - r) * cr2);
   Int_t i = TMath::LocMin(2, saf);
   if (TGeoShape::IsCloseToPhi(saf[i], point, c1, s1, c2, s2)) {
      TGeoShape::NormalPhi(point, dir, norm, c1, s1, c2, s2);
      return;
   }

   Double_t phi = TMath::ATan2(point[1], point[0]);
   Double_t cphi = TMath::Cos(phi);
   Double_t sphi = TMath::Sin(phi);

   if (i == 0) {
      norm[0] = cr1 * cphi;
      norm[1] = cr1 * sphi;
      norm[2] = -tg1 * cr1;
   } else {
      norm[0] = cr2 * cphi;
      norm[1] = cr2 * sphi;
      norm[2] = -tg2 * cr2;
   }

   if (norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2] < 0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
      norm[2] = -norm[2];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute distance to a conical surface with :
/// - r1, z1 - radius and Z position of lower base
/// - r2, z2 - radius and Z position of upper base
/// - phi1, phi2 - phi limits

Double_t TGeoVGConeSeg::DistToCons(const Double_t *point, const Double_t *dir, Double_t r1, Double_t z1, Double_t r2,
                                 Double_t z2, Double_t phi1, Double_t phi2)
{
   Double_t dz = z2 - z1;
   if (dz <= 0) {
      return TGeoShape::Big();
   }

   Double_t dphi = phi2 - phi1;
   Bool_t hasphi = kTRUE;
   if (dphi >= 360.)
      hasphi = kFALSE;
   if (dphi < 0)
      dphi += 360.;
   //   printf("phi1=%f phi2=%f dphi=%f\n", phi1, phi2, dphi);

   Double_t ro0 = 0.5 * (r1 + r2);
   Double_t fz = (r2 - r1) / dz;
   Double_t r0sq = point[0] * point[0] + point[1] * point[1];
   Double_t rc = ro0 + fz * (point[2] - 0.5 * (z1 + z2));

   Double_t a = dir[0] * dir[0] + dir[1] * dir[1] - fz * fz * dir[2] * dir[2];
   Double_t b = point[0] * dir[0] + point[1] * dir[1] - fz * rc * dir[2];
   Double_t c = r0sq - rc * rc;

   if (a == 0)
      return TGeoShape::Big();
   a = 1. / a;
   b *= a;
   c *= a;
   Double_t delta = b * b - c;
   if (delta < 0)
      return TGeoShape::Big();
   delta = TMath::Sqrt(delta);

   Double_t snxt = -b - delta;
   Double_t ptnew[3];
   Double_t ddp, phi;
   if (snxt > 0) {
      // check Z range
      ptnew[2] = point[2] + snxt * dir[2];
      if (((ptnew[2] - z1) * (ptnew[2] - z2)) < 0) {
         // check phi range
         if (!hasphi)
            return snxt;
         ptnew[0] = point[0] + snxt * dir[0];
         ptnew[1] = point[1] + snxt * dir[1];
         phi = TMath::ATan2(ptnew[1], ptnew[0]) * TMath::RadToDeg();
         if (phi < 0)
            phi += 360.;
         ddp = phi - phi1;
         if (ddp < 0)
            ddp += 360.;
         // printf("snxt1=%f phi=%f ddp=%f\n", snxt, phi, ddp);
         if (ddp <= dphi)
            return snxt;
      }
   }
   snxt = -b + delta;
   if (snxt > 0) {
      // check Z range
      ptnew[2] = point[2] + snxt * dir[2];
      if (((ptnew[2] - z1) * (ptnew[2] - z2)) < 0) {
         // check phi range
         if (!hasphi)
            return snxt;
         ptnew[0] = point[0] + snxt * dir[0];
         ptnew[1] = point[1] + snxt * dir[1];
         phi = TMath::ATan2(ptnew[1], ptnew[0]) * TMath::RadToDeg();
         if (phi < 0)
            phi += 360.;
         ddp = phi - phi1;
         if (ddp < 0)
            ddp += 360.;
         // printf("snxt2=%f phi=%f ddp=%f\n", snxt, phi, ddp);
         if (ddp <= dphi)
            return snxt;
      }
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Static: compute distance from inside point to surface of the tube segment

Double_t TGeoVGConeSeg::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                      Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1, Double_t s1,
                                      Double_t c2, Double_t s2, Double_t cm, Double_t sm, Double_t cdfi)
{
   if (dz <= 0)
      return TGeoShape::Big();
   // Do Z
   Double_t scone = TGeoVGCone::DistFromInsideS(point, dir, dz, rmin1, rmax1, rmin2, rmax2);
   if (scone <= 0)
      return 0.0;
   Double_t sfmin = TGeoShape::Big();
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t cpsi = point[0] * cm + point[1] * sm;
   if (cpsi > r * cdfi + TGeoShape::Tolerance()) {
      sfmin = TGeoShape::DistToPhiMin(point, dir, s1, c1, s2, c2, sm, cm);
      return TMath::Min(scone, sfmin);
   }
   // Point on the phi boundary or outside
   // which one: phi1 or phi2
   Double_t ddotn, xi, yi;
   if (TMath::Abs(point[1] - s1 * r) < TMath::Abs(point[1] - s2 * r)) {
      ddotn = s1 * dir[0] - c1 * dir[1];
      if (ddotn >= 0)
         return 0.0;
      ddotn = -s2 * dir[0] + c2 * dir[1];
      if (ddotn <= 0)
         return scone;
      sfmin = s2 * point[0] - c2 * point[1];
      if (sfmin <= 0)
         return scone;
      sfmin /= ddotn;
      if (sfmin >= scone)
         return scone;
      xi = point[0] + sfmin * dir[0];
      yi = point[1] + sfmin * dir[1];
      if (yi * cm - xi * sm < 0)
         return scone;
      return sfmin;
   }
   ddotn = -s2 * dir[0] + c2 * dir[1];
   if (ddotn >= 0)
      return 0.0;
   ddotn = s1 * dir[0] - c1 * dir[1];
   if (ddotn <= 0)
      return scone;
   sfmin = -s1 * point[0] + c1 * point[1];
   if (sfmin <= 0)
      return scone;
   sfmin /= ddotn;
   if (sfmin >= scone)
      return scone;
   xi = point[0] + sfmin * dir[0];
   yi = point[1] + sfmin * dir[1];
   if (yi * cm - xi * sm > 0)
      return scone;
   return sfmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: compute distance from outside point to surface of arbitrary tube

Double_t TGeoVGConeSeg::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                       Double_t rmax1, Double_t rmin2, Double_t rmax2, Double_t c1, Double_t s1,
                                       Double_t c2, Double_t s2, Double_t cm, Double_t sm, Double_t cdfi)
{
   if (dz <= 0)
      return TGeoShape::Big();
   Double_t r2, cpsi;
   // check Z planes
   Double_t xi, yi, zi;
   Double_t b, delta;
   zi = dz - TMath::Abs(point[2]);
   Double_t rin, rout;
   Double_t s = TGeoShape::Big();
   Double_t snxt = TGeoShape::Big();
   Bool_t in = kFALSE;
   Bool_t inz = (zi < 0) ? kFALSE : kTRUE;
   if (!inz) {
      if (point[2] * dir[2] >= 0)
         return TGeoShape::Big();
      s = -zi / TMath::Abs(dir[2]);
      xi = point[0] + s * dir[0];
      yi = point[1] + s * dir[1];
      r2 = xi * xi + yi * yi;
      if (dir[2] > 0) {
         rin = rmin1;
         rout = rmax1;
      } else {
         rin = rmin2;
         rout = rmax2;
      }
      if ((rin * rin <= r2) && (r2 <= rout * rout)) {
         cpsi = xi * cm + yi * sm;
         if (cpsi >= (cdfi * TMath::Sqrt(r2)))
            return s;
      }
   }
   Double_t zinv = 1. / dz;
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t ro1 = 0.5 * (rmin1 + rmin2);
   Bool_t hasrmin = (ro1 > 0) ? kTRUE : kFALSE;
   Double_t tg1 = 0.0;
   Bool_t inrmin = kFALSE;
   rin = 0.0;
   if (hasrmin) {
      tg1 = 0.5 * (rmin2 - rmin1) * zinv;
      rin = ro1 + tg1 * point[2];
      if (rsq > rin * (rin - TGeoShape::Tolerance()))
         inrmin = kTRUE;
   } else {
      inrmin = kTRUE;
   }
   Double_t ro2 = 0.5 * (rmax1 + rmax2);
   Double_t tg2 = 0.5 * (rmax2 - rmax1) * zinv;
   rout = ro2 + tg2 * point[2];
   Bool_t inrmax = kFALSE;
   if (r < rout + TGeoShape::Tolerance())
      inrmax = kTRUE;
   Bool_t inphi = kFALSE;
   cpsi = point[0] * cm + point[1] * sm;
   if (cpsi > r * cdfi - TGeoShape::Tolerance())
      inphi = kTRUE;
   in = inz & inrmin & inrmax & inphi;
   // If inside, we are most likely on a boundary within machine precision.
   if (in) {
      Double_t safphi = (cpsi - r * cdfi) * TMath::Sqrt(1. - cdfi * cdfi);
      Double_t safrmin = (hasrmin) ? TMath::Abs(r - rin) : (TGeoShape::Big());
      Double_t safrmax = TMath::Abs(r - rout);
      // check if on Z boundaries
      if (zi < safrmax && zi < safrmin && zi < safphi) {
         if (point[2] * dir[2] < 0)
            return 0.0;
         return TGeoShape::Big();
      }
      // check if on Rmax boundary
      if (safrmax < safrmin && safrmax < safphi) {
         Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] - tg2 * dir[2] * r;
         if (ddotn <= 0)
            return 0.0;
         return TGeoShape::Big();
      }
      // check if on phi boundary
      if (safphi < safrmin) {
         // We may cross again a phi of rmin boundary
         // check first if we are on phi1 or phi2
         Double_t un;
         if (point[0] * c1 + point[1] * s1 > point[0] * c2 + point[1] * s2) {
            un = dir[0] * s1 - dir[1] * c1;
            if (un < 0)
               return 0.0;
            if (cdfi >= 0)
               return TGeoShape::Big();
            un = -dir[0] * s2 + dir[1] * c2;
            if (un < 0) {
               s = -point[0] * s2 + point[1] * c2;
               if (s > 0) {
                  s /= (-un);
                  zi = point[2] + s * dir[2];
                  if (TMath::Abs(zi) <= dz) {
                     xi = point[0] + s * dir[0];
                     yi = point[1] + s * dir[1];
                     if ((yi * cm - xi * sm) > 0) {
                        r2 = xi * xi + yi * yi;
                        rin = ro1 + tg1 * zi;
                        rout = ro2 + tg2 * zi;
                        if ((rin * rin <= r2) && (rout * rout >= r2))
                           return s;
                     }
                  }
               }
            }
         } else {
            un = -dir[0] * s2 + dir[1] * c2;
            if (un < 0)
               return 0.0;
            if (cdfi >= 0)
               return TGeoShape::Big();
            un = dir[0] * s1 - dir[1] * c1;
            if (un < 0) {
               s = point[0] * s1 - point[1] * c1;
               if (s > 0) {
                  s /= (-un);
                  zi = point[2] + s * dir[2];
                  if (TMath::Abs(zi) <= dz) {
                     xi = point[0] + s * dir[0];
                     yi = point[1] + s * dir[1];
                     if ((yi * cm - xi * sm) < 0) {
                        r2 = xi * xi + yi * yi;
                        rin = ro1 + tg1 * zi;
                        rout = ro2 + tg2 * zi;
                        if ((rin * rin <= r2) && (rout * rout >= r2))
                           return s;
                     }
                  }
               }
            }
         }
         // We may also cross rmin, second solution coming from outside
         Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] - tg1 * dir[2] * r;
         if (ddotn >= 0)
            return TGeoShape::Big();
         if (cdfi >= 0)
            return TGeoShape::Big();
         TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
         if (delta < 0)
            return TGeoShape::Big();
         snxt = -b - delta;
         if (snxt < 0)
            return TGeoShape::Big();
         snxt = -b + delta;
         zi = point[2] + snxt * dir[2];
         if (TMath::Abs(zi) > dz)
            return TGeoShape::Big();
         xi = point[0] + snxt * dir[0];
         yi = point[1] + snxt * dir[1];
         r2 = xi * xi + yi * yi;
         cpsi = xi * cm + yi * sm;
         if (cpsi >= (cdfi * TMath::Sqrt(r2)))
            return snxt;
         return TGeoShape::Big();
      }
      // We are on rmin boundary: we may cross again rmin or a phi facette
      Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] - tg1 * dir[2] * r;
      if (ddotn >= 0)
         return 0.0;
      TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
      if (delta < 0)
         return 0.0;
      snxt = -b + delta;
      if (snxt < 0)
         return TGeoShape::Big();
      if (TMath::Abs(-b - delta) > snxt)
         return TGeoShape::Big();
      zi = point[2] + snxt * dir[2];
      if (TMath::Abs(zi) > dz)
         return TGeoShape::Big();
      // OK, we cross rmin at snxt - check if within phi range
      xi = point[0] + snxt * dir[0];
      yi = point[1] + snxt * dir[1];
      r2 = xi * xi + yi * yi;
      cpsi = xi * cm + yi * sm;
      if (cpsi >= (cdfi * TMath::Sqrt(r2)))
         return snxt;
      // we cross rmin in the phi gap - we may cross a phi facette
      if (cdfi >= 0)
         return TGeoShape::Big();
      Double_t un = -dir[0] * s1 + dir[1] * c1;
      if (un > 0) {
         s = point[0] * s1 - point[1] * c1;
         if (s >= 0) {
            s /= un;
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               if ((yi * cm - xi * sm) <= 0) {
                  r2 = xi * xi + yi * yi;
                  rin = ro1 + tg1 * zi;
                  rout = ro2 + tg2 * zi;
                  if ((rin * rin <= r2) && (rout * rout >= r2))
                     return s;
               }
            }
         }
      }
      un = dir[0] * s2 - dir[1] * c2;
      if (un > 0) {
         s = (point[1] * c2 - point[0] * s2) / un;
         if (s >= 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               if ((yi * cm - xi * sm) >= 0) {
                  r2 = xi * xi + yi * yi;
                  rin = ro1 + tg1 * zi;
                  rout = ro2 + tg2 * zi;
                  if ((rin * rin <= r2) && (rout * rout >= r2))
                     return s;
               }
            }
         }
      }
      return TGeoShape::Big();
   }

   // The point is really outside
   Double_t sr1 = TGeoShape::Big();
   if (!inrmax) {
      // check crossing with outer cone
      TGeoVGCone::DistToCone(point, dir, dz, rmax1, rmax2, b, delta);
      if (delta >= 0) {
         s = -b - delta;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               r2 = xi * xi + yi * yi;
               cpsi = xi * cm + yi * sm;
               if (cpsi >= (cdfi * TMath::Sqrt(r2)))
                  return s; // rmax crossing
            }
         }
         s = -b + delta;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               r2 = xi * xi + yi * yi;
               cpsi = xi * cm + yi * sm;
               if (cpsi >= (cdfi * TMath::Sqrt(r2)))
                  sr1 = s;
            }
         }
      }
   }
   // check crossing with inner cone
   Double_t sr2 = TGeoShape::Big();
   TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
   if (delta >= 0) {
      s = -b - delta;
      if (s > 0) {
         zi = point[2] + s * dir[2];
         if (TMath::Abs(zi) <= dz) {
            xi = point[0] + s * dir[0];
            yi = point[1] + s * dir[1];
            r2 = xi * xi + yi * yi;
            cpsi = xi * cm + yi * sm;
            if (cpsi >= (cdfi * TMath::Sqrt(r2)))
               sr2 = s;
         }
      }
      if (sr2 > 1E10) {
         s = -b + delta;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               r2 = xi * xi + yi * yi;
               cpsi = xi * cm + yi * sm;
               if (cpsi >= (cdfi * TMath::Sqrt(r2)))
                  sr2 = s;
            }
         }
      }
   }
   snxt = TMath::Min(sr1, sr2);
   // Check phi crossing
   s = TGeoShape::DistToPhiMin(point, dir, s1, c1, s2, c2, sm, cm, kFALSE);
   if (s > snxt)
      return snxt;
   zi = point[2] + s * dir[2];
   if (TMath::Abs(zi) > dz)
      return snxt;
   xi = point[0] + s * dir[0];
   yi = point[1] + s * dir[1];
   r2 = xi * xi + yi * yi;
   rout = ro2 + tg2 * zi;
   if (r2 > rout * rout)
      return snxt;
   rin = ro1 + tg1 * zi;
   if (r2 >= rin * rin)
      return s; // phi crossing
   return snxt;
}

// ////////////////////////////////////////////////////////////////////////////////
// /// compute closest distance from point px,py to each corner

// Int_t TGeoVGConeSeg::DistancetoPrimitive(Int_t px, Int_t py)
// {
//    Int_t n = gGeoManager->GetNsegments() + 1;
//    const Int_t numPoints = 4 * n;
//    return ShapeDistancetoPrimitive(numPoints, px, py);
// }

////////////////////////////////////////////////////////////////////////////////
/// Divide this cone segment shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. For Z division
/// creates all volumes with different shapes and returns pointer to volume that
/// was divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *
TGeoVGConeSeg::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Double_t dphi;
   Int_t id;
   Double_t end = start + ndiv * step;
   switch (iaxis) {
   case 1: //---               R division
      Error("Divide", "division of a cone segment on R not implemented");
      return nullptr;
   case 2: //---               Phi division
      dphi = GetDPhi() * TMath::RadToDeg();
      if (dphi < 0)
         dphi += 360.;
      finder = new TGeoPatternCylPhi(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGConeSeg(GetDz(), GetRmin1(), GetRmax1(), GetRmin2(), GetRmax2(), -step / 2, step / 2);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      opt = "Phi";
      for (id = 0; id < ndiv; id++) {
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 3: //---                 Z division
      finder = new TGeoPatternZ(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         Double_t z1 = start + id * step;
         Double_t z2 = start + (id + 1) * step;
         Double_t rmin1n = 0.5 * (GetRmin1() * (GetDz() - z1) + GetRmin2() * (GetDz() + z1)) / GetDz();
         Double_t rmax1n = 0.5 * (GetRmax1() * (GetDz() - z1) + GetRmax2() * (GetDz() + z1)) / GetDz();
         Double_t rmin2n = 0.5 * (GetRmin1() * (GetDz() - z2) + GetRmin2() * (GetDz() + z2)) / GetDz();
         Double_t rmax2n = 0.5 * (GetRmax1() * (GetDz() - z2) + GetRmax2() * (GetDz() + z2)) / GetDz();
         shape = new TGeoVGConeSeg(step / 2, rmin1n, rmax1n, rmin2n, rmax2n, GetSPhi() * TMath::RadToDeg(), (GetSPhi() + GetDPhi()) * TMath::RadToDeg());
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Z";
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   default: Error("Divide", "Wrong axis type for division"); return nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGConeSeg::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 2:
      xlo = GetSPhi() * TMath::RadToDeg();
      xhi = (GetSPhi() + GetDPhi()) * TMath::RadToDeg();
      dx = xhi - xlo;
      return dx;
   case 3:
      xlo = -GetDz();
      xhi = GetDz();
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGConeSeg::GetBoundingCylinder(Double_t *param) const
{
   param[0] = TMath::Min(GetRmin1(), GetRmin2()); // Rmin
   param[0] *= param[0];
   param[1] = TMath::Max(GetRmax1(), GetRmax2()); // Rmax
   param[1] *= param[1];
   param[2] = ((GetSPhi() * TMath::RadToDeg()) < 0) ? ((GetSPhi() * TMath::RadToDeg()) + 360.) : GetSPhi() * TMath::RadToDeg(); // Phi1
   param[3] = (GetSPhi() + GetDPhi()) * TMath::RadToDeg();  // Phi2
   while (param[3] < param[2])
      param[3] += 360.;
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGConeSeg::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (!mother->TestShapeBit(kGeoConeSeg)) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t rmin1, rmax1, rmin2, rmax2, dz;
   rmin1 = GetRmin1();
   rmax1 = GetRmax1();
   rmin2 = GetRmin2();
   rmax2 = GetRmax2();
   dz = GetDz();
   if (GetDz() < 0)
      dz = ((TGeoVGConeSeg *)mother)->GetDz();
   if (GetRmin1() < 0)
      rmin1 = ((TGeoVGConeSeg *)mother)->GetRmin1();
   if ((GetRmax1() < 0) || (GetRmax1() < GetRmin1()))
      rmax1 = ((TGeoVGConeSeg *)mother)->GetRmax1();
   if (GetRmin2() < 0)
      rmin2 = ((TGeoVGConeSeg *)mother)->GetRmin2();
   if ((GetRmax2() < 0) || (GetRmax2() < GetRmin2()))
      rmax2 = ((TGeoVGConeSeg *)mother)->GetRmax2();

   return (new TGeoVGConeSeg(GetName(), dz, rmin1, rmax1, rmin2, rmax2, GetSPhi() * TMath::RadToDeg(), (GetSPhi() + GetDPhi()) * TMath::RadToDeg()));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGConeSeg::InspectShape() const
{
   printf("*** Shape %s: TGeoVGConeSeg ***\n", GetName());
   printf("    dz    = %11.5f\n", GetDz());
   printf("    Rmin1 = %11.5f\n", GetRmin1());
   printf("    Rmax1 = %11.5f\n", GetRmax1());
   printf("    Rmin2 = %11.5f\n", GetRmin2());
   printf("    Rmax2 = %11.5f\n", GetRmax2());
   printf("    phi1  = %11.5f\n", GetSPhi() * TMath::RadToDeg());
   printf("    phi2  = %11.5f\n", (GetSPhi() + GetDPhi()) * TMath::RadToDeg());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

///////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGConeSeg::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t nbPnts = 4 * n;
   Int_t nbSegs = 2 * nbPnts;
   Int_t nbPols = nbPnts - 2;
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

void TGeoVGConeSeg::SetSegsAndPols(TBuffer3D &buffer) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t c = GetBasicColor();

   memset(buffer.fSegs, 0, buffer.NbSegs() * 3 * sizeof(Int_t));
   for (i = 0; i < 4; i++) {
      for (j = 1; j < n; j++) {
         buffer.fSegs[(i * n + j - 1) * 3] = c;
         buffer.fSegs[(i * n + j - 1) * 3 + 1] = i * n + j - 1;
         buffer.fSegs[(i * n + j - 1) * 3 + 2] = i * n + j;
      }
   }
   for (i = 4; i < 6; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i * n + j) * 3] = c + 1;
         buffer.fSegs[(i * n + j) * 3 + 1] = (i - 4) * n + j;
         buffer.fSegs[(i * n + j) * 3 + 2] = (i - 2) * n + j;
      }
   }
   for (i = 6; i < 8; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i * n + j) * 3] = c;
         buffer.fSegs[(i * n + j) * 3 + 1] = 2 * (i - 6) * n + j;
         buffer.fSegs[(i * n + j) * 3 + 2] = (2 * (i - 6) + 1) * n + j;
      }
   }

   Int_t indx = 0;
   memset(buffer.fPols, 0, buffer.NbPols() * 6 * sizeof(Int_t));
   i = 0;
   for (j = 0; j < n - 1; j++) {
      buffer.fPols[indx++] = c;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = (4 + i) * n + j + 1;
      buffer.fPols[indx++] = (2 + i) * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j;
      buffer.fPols[indx++] = i * n + j;
   }
   i = 1;
   for (j = 0; j < n - 1; j++) {
      buffer.fPols[indx++] = c;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = i * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j;
      buffer.fPols[indx++] = (2 + i) * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j + 1;
   }
   i = 2;
   for (j = 0; j < n - 1; j++) {
      buffer.fPols[indx++] = c + i;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = (i - 2) * 2 * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j;
      buffer.fPols[indx++] = ((i - 2) * 2 + 1) * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j + 1;
   }
   i = 3;
   for (j = 0; j < n - 1; j++) {
      buffer.fPols[indx++] = c + i;
      buffer.fPols[indx++] = 4;
      buffer.fPols[indx++] = (4 + i) * n + j + 1;
      buffer.fPols[indx++] = ((i - 2) * 2 + 1) * n + j;
      buffer.fPols[indx++] = (4 + i) * n + j;
      buffer.fPols[indx++] = (i - 2) * 2 * n + j;
   }
   buffer.fPols[indx++] = c + 2;
   buffer.fPols[indx++] = 4;
   buffer.fPols[indx++] = 6 * n;
   buffer.fPols[indx++] = 4 * n;
   buffer.fPols[indx++] = 7 * n;
   buffer.fPols[indx++] = 5 * n;
   buffer.fPols[indx++] = c + 2;
   buffer.fPols[indx++] = 4;
   buffer.fPols[indx++] = 6 * n - 1;
   buffer.fPols[indx++] = 8 * n - 1;
   buffer.fPols[indx++] = 5 * n - 1;
   buffer.fPols[indx++] = 7 * n - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute the closest distance from given point to this shape.

Double_t TGeoVGConeSeg::SafetyS(const Double_t *point, Bool_t in, Double_t dz, Double_t rmin1, Double_t rmax1,
                              Double_t rmin2, Double_t rmax2, Double_t phi1, Double_t phi2, Int_t skipz)
{
   Double_t safe = TGeoVGCone::SafetyS(point, in, dz, rmin1, rmax1, rmin2, rmax2, skipz);
   if ((phi2 - phi1) >= 360.)
      return safe;
   Double_t safphi = TGeoShape::SafetyPhi(point, in, phi1, phi2);
   if (in)
      return TMath::Min(safe, safphi);
   if (safe > 1.E10)
      return safphi;
   return TMath::Max(safe, safphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGConeSeg::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz    = " << GetDz() << ";" << std::endl;
   out << "   rmin1 = " << GetRmin1() << ";" << std::endl;
   out << "   rmax1 = " << GetRmax1() << ";" << std::endl;
   out << "   rmin2 = " << GetRmin2() << ";" << std::endl;
   out << "   rmax2 = " << GetRmax2() << ";" << std::endl;
   out << "   phi1  = " << GetSPhi() * TMath::RadToDeg() << ";" << std::endl;
   out << "   phi2  = " << (GetSPhi() + GetDPhi()) * TMath::RadToDeg() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGConeSeg(\"" << GetName()
       << "\", dz,rmin1,rmax1,rmin2,rmax2,phi1,phi2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the cone segment.

void TGeoVGConeSeg::SetConsDimensions(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2,
                                    Double_t phi1, Double_t phi2)
{
   SetDz(dz);
   SetRmin1(rmin1);
   SetRmax1(rmax1);
   SetRmin2(rmin2);
   SetRmax2(rmax2);
   auto newPhi1 = phi1;
   while (newPhi1 < 0)
      newPhi1 += 360.;
   SetSPhi(newPhi1 * TMath::DegToRad());
   auto newPhi2 = phi2;
   while (newPhi2 <= newPhi1)
      newPhi2 += 360.;
   SetDPhi((newPhi2 - newPhi1) * TMath::DegToRad());

   if (TGeoShape::IsSameWithinTolerance(GetSPhi()* TMath::RadToDeg(), (GetSPhi() + GetDPhi()) * TMath::RadToDeg()))
      Fatal("SetConsDimensions", "In shape %s invalid phi1=%g, phi2=%g\n", GetName(), phi1, phi2);
   InitTrigonometry();
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the cone segment from an array.

void TGeoVGConeSeg::SetDimensions(Double_t *param)
{
   Double_t dz = param[0];
   Double_t rmin1 = param[1];
   Double_t rmax1 = param[2];
   Double_t rmin2 = param[3];
   Double_t rmax2 = param[4];
   Double_t phi1 = param[5];
   Double_t phi2 = param[6];
   SetConsDimensions(dz, rmin1, rmax1, rmin2, rmax2, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Create cone segment mesh points.

void TGeoVGConeSeg::SetPoints(Double_t *points) const
{
   Int_t j, n;
   Float_t dphi, phi, phi1, phi2, dz;

   n = gGeoManager->GetNsegments() + 1;
   dz = GetDz();
   phi1 = GetSPhi() * TMath::RadToDeg();
   phi2 = (GetSPhi() + GetDPhi()) * TMath::RadToDeg();

   dphi = (phi2 - phi1) / (n - 1);

   Int_t indx = 0;

   if (points) {
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmin1() * TMath::Cos(phi);
         points[indx++] = GetRmin1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmax1() * TMath::Cos(phi);
         points[indx++] = GetRmax1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmin2() * TMath::Cos(phi);
         points[indx++] = GetRmin2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmax2() * TMath::Cos(phi);
         points[indx++] = GetRmax2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create cone segment mesh points.

void TGeoVGConeSeg::SetPoints(Float_t *points) const
{
   Int_t j, n;
   Float_t dphi, phi, phi1, phi2, dz;

   n = gGeoManager->GetNsegments() + 1;
   dz = GetDz();
   phi1 = GetSPhi() * TMath::RadToDeg();
   phi2 = (GetSPhi() + GetDPhi()) * TMath::RadToDeg();

   dphi = (phi2 - phi1) / (n - 1);

   Int_t indx = 0;

   if (points) {
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmin1() * TMath::Cos(phi);
         points[indx++] = GetRmin1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmax1() * TMath::Cos(phi);
         points[indx++] = GetRmax1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmin2() * TMath::Cos(phi);
         points[indx++] = GetRmin2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
      for (j = 0; j < n; j++) {
         phi = (GetSPhi() * TMath::RadToDeg() + j * dphi) * TMath::DegToRad();
         points[indx++] = GetRmax2() * TMath::Cos(phi);
         points[indx++] = GetRmax2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGConeSeg::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   nvert = n * 4;
   nsegs = n * 8;
   npols = n * 4 - 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGConeSeg::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t numPoints = n * 4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill size of this 3-D object

void TGeoVGConeSeg::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGConeSeg::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

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

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoVGConeSeg::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints / 2) * 2) {
      Error("GetPointsOnSegments", "Npoints must be even number");
      return kFALSE;
   }
   Int_t nc = (Int_t)TMath::Sqrt(0.5 * npoints);
   Double_t dphi = GetDPhi() / (nc - 1);
   Double_t phi = 0;
   Double_t phi1 = GetSPhi();
   Int_t ntop = npoints / 2 - nc * (nc - 1);
   Double_t dz = 2 * GetDz() / (nc - 1);
   Double_t z = 0;
   Double_t rmin = 0.;
   Double_t rmax = 0.;
   Int_t icrt = 0;
   Int_t nphi = nc;
   // loop z sections
   for (Int_t i = 0; i < nc; i++) {
      if (i == (nc - 1)) {
         nphi = ntop;
         dphi = GetDPhi() / (nphi - 1);
      }
      z = -GetDz() + i * dz;
      rmin = 0.5 * (GetRmin1() + GetRmin2()) + 0.5 * (GetRmin2() - GetRmin1()) * z / GetDz();
      rmax = 0.5 * (GetRmax1() + GetRmax2()) + 0.5 * (GetRmax2() - GetRmax1()) * z / GetDz();
      // loop points on circle sections
      for (Int_t j = 0; j < nphi; j++) {
         phi = phi1 + j * dphi;
         array[icrt++] = rmin * TMath::Cos(phi);
         array[icrt++] = rmin * TMath::Sin(phi);
         array[icrt++] = z;
         array[icrt++] = rmax * TMath::Cos(phi);
         array[icrt++] = rmax * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

#endif
