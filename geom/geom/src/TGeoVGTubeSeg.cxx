// @(#)root/geom:$Id$

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoVGTubeSeg.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGTubeSeg);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGTubeSeg::TGeoVGTubeSeg()
   : Base_t("", 0., 0., 0., 0., 0.),
     fS1(0.), fC1(0.), fS2(0.), fC2(0.), fSm(0.), fCm(0.), fCdfi(0.)
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius.
/// The segment will be from phiStart to phiEnd expressed in degree.

TGeoVGTubeSeg::TGeoVGTubeSeg(Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart, Double_t phiEnd)
   : Base_t("", rmin, rmax, dz, phiStart * TMath::DegToRad(), (phiEnd - phiStart) * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetTubsDimensions(rmin, rmax, dz, phiStart, phiEnd);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
/// The segment will be from phiStart to phiEnd expressed in degree.

TGeoVGTubeSeg::TGeoVGTubeSeg(const char *name, Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart,
                         Double_t phiEnd)
   : Base_t(name, rmin, rmax, dz, phiStart * TMath::DegToRad(), (phiEnd - phiStart) * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetTubsDimensions(rmin, rmax, dz, phiStart, phiEnd);
   SetName(name);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = Rmin
///  - param[1] = Rmax
///  - param[2] = dz
///  - param[3] = phi1
///  - param[4] = phi2

TGeoVGTubeSeg::TGeoVGTubeSeg(Double_t *param) 
  : Base_t("", param[0], param[1], param[2], param[3], param[4] - param[3])
{
   SetShapeBit(TGeoShape::kGeoTubeSeg);
   SetDimensions(param);
   ComputeBBox();
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Init frequently used trigonometric values

void TGeoVGTubeSeg::InitTrigonometry()
{
   Double_t phi1 = sphi();
   Double_t phi2 = sphi() + dphi();
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
/// compute bounding box of the tube segment

void TGeoVGTubeSeg::ComputeBBox()
{
   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax() * fC1;
   yc[0] = rmax() * fS1;
   xc[1] = rmax() * fC2;
   yc[1] = rmax() * fS2;
   xc[2] = rmin() * fC1;
   yc[2] = rmin() * fS1;
   xc[3] = rmin() * fC2;
   yc[3] = rmin() * fS2;

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];

   Double_t dp = dphi() * TMath::RadToDeg();
   if (dp < 0)
      dp += 360;
   Double_t ddp = -sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      xmax = rmax();
   ddp = 90 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      ymax = rmax();
   ddp = 180 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      xmin = -rmax();
   ddp = 270 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      ymin = -rmax();
   fOrigin[0] = (xmax + xmin) / 2;
   fOrigin[1] = (ymax + ymin) / 2;
   fOrigin[2] = 0;
   fDX = (xmax - xmin) / 2;
   fDY = (ymax - ymin) / 2;
   fDZ = z();
}


////////////////////////////////////////////////////////////////////////////////
/// Static: Computes capacity of the shape in [length^3]

Double_t TGeoVGTubeSeg::Capacity(Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart, Double_t phiEnd)
{
   Double_t capacity = TMath::Abs(phiEnd - phiStart) * TMath::DegToRad() * (rmax * rmax - rmin * rmin) * dz;
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute normal to closest surface from POINT.

void TGeoVGTubeSeg::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t rmin,
                                 Double_t rmax, Double_t /*dz*/, Double_t c1, Double_t s1, Double_t c2, Double_t s2)
{
   Double_t saf[2];
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);
   saf[0] = (rmin > 1E-10) ? TMath::Abs(r - rmin) : TGeoShape::Big();
   saf[1] = TMath::Abs(rmax - r);
   Int_t i = TMath::LocMin(2, saf);
   if (TGeoShape::IsCloseToPhi(saf[i], point, c1, s1, c2, s2)) {
      TGeoShape::NormalPhi(point, dir, norm, c1, s1, c2, s2);
      return;
   }
   norm[2] = 0;
   Double_t phi = TMath::ATan2(point[1], point[0]);
   norm[0] = TMath::Cos(phi);
   norm[1] = TMath::Sin(phi);
   if (norm[0] * dir[0] + norm[1] * dir[1] < 0) {
      norm[0] = -norm[0];
      norm[1] = -norm[1];
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute distance from inside point to surface of the tube segment (static)
/// Boundary safe algorithm.
/// Do Z

Double_t TGeoVGTubeSeg::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax,
                                      Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2, Double_t cm,
                                      Double_t sm, Double_t cdfi)
{
   Double_t stube = TGeoTube::DistFromInsideS(point, dir, rmin, rmax, dz);
   if (stube <= 0)
      return 0.0;
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t cpsi = point[0] * cm + point[1] * sm;
   if (cpsi > r * cdfi + TGeoShape::Tolerance()) {
      Double_t sfmin = TGeoShape::DistToPhiMin(point, dir, s1, c1, s2, c2, sm, cm);
      return TMath::Min(stube, sfmin);
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
         return stube;
      Double_t sfmin = s2 * point[0] - c2 * point[1];
      if (sfmin <= 0)
         return stube;
      sfmin /= ddotn;
      if (sfmin >= stube)
         return stube;
      xi = point[0] + sfmin * dir[0];
      yi = point[1] + sfmin * dir[1];
      if (yi * cm - xi * sm < 0)
         return stube;
      return sfmin;
   }
   ddotn = -s2 * dir[0] + c2 * dir[1];
   if (ddotn >= 0)
      return 0.0;
   ddotn = s1 * dir[0] - c1 * dir[1];
   if (ddotn <= 0)
      return stube;
   Double_t sfmin = -s1 * point[0] + c1 * point[1];
   if (sfmin <= 0)
      return stube;
   sfmin /= ddotn;
   if (sfmin >= stube)
      return stube;
   xi = point[0] + sfmin * dir[0];
   yi = point[1] + sfmin * dir[1];
   if (yi * cm - xi * sm > 0)
      return stube;
   return sfmin;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Static method to compute distance to arbitrary tube segment from outside point
/// Boundary safe algorithm.

Double_t TGeoVGTubeSeg::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax,
                                       Double_t dz, Double_t c1, Double_t s1, Double_t c2, Double_t s2, Double_t cm,
                                       Double_t sm, Double_t cdfi)
{
   Double_t r2, cpsi, s;
   // check Z planes
   Double_t xi, yi, zi;
   zi = dz - TMath::Abs(point[2]);
   Double_t rmaxsq = rmax * rmax;
   Double_t rminsq = rmin * rmin;
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
      if ((rminsq <= r2) && (r2 <= rmaxsq)) {
         cpsi = (xi * cm + yi * sm) / TMath::Sqrt(r2);
         if (cpsi >= cdfi)
            return s;
      }
   }

   // check outer cyl. surface
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);
   Double_t nsq = dir[0] * dir[0] + dir[1] * dir[1];
   Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
   Double_t b, d;
   Bool_t inrmax = kFALSE;
   Bool_t inrmin = kFALSE;
   Bool_t inphi = kFALSE;
   if (rsq <= rmaxsq + TGeoShape::Tolerance())
      inrmax = kTRUE;
   if (rsq >= rminsq - TGeoShape::Tolerance())
      inrmin = kTRUE;
   cpsi = point[0] * cm + point[1] * sm;
   if (cpsi > r * cdfi - TGeoShape::Tolerance())
      inphi = kTRUE;
   in = inz & inrmin & inrmax & inphi;
   // If inside, we are most likely on a boundary within machine precision.
   if (in) {
      Bool_t checkout = kFALSE;
      Double_t safphi = (cpsi - r * cdfi) * TMath::Sqrt(1. - cdfi * cdfi);
      //      Double_t sch, cch;
      // check if on Z boundaries
      if (zi < rmax - r) {
         if (TGeoShape::IsSameWithinTolerance(rmin, 0) || (zi < r - rmin)) {
            if (zi < safphi) {
               if (point[2] * dir[2] < 0)
                  return 0.0;
               return TGeoShape::Big();
            }
         }
      }
      if ((rmaxsq - rsq) < (rsq - rminsq))
         checkout = kTRUE;
      // check if on Rmax boundary
      if (checkout && (rmax - r < safphi)) {
         if (rdotn >= 0)
            return TGeoShape::Big();
         return 0.0;
      }
      if (TMath::Abs(nsq) < TGeoShape::Tolerance())
         return TGeoShape::Big();
      // check if on phi boundary
      if (TGeoShape::IsSameWithinTolerance(rmin, 0) || (safphi < r - rmin)) {
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
                     r2 = xi * xi + yi * yi;
                     if ((rminsq <= r2) && (r2 <= rmaxsq)) {
                        if ((yi * cm - xi * sm) > 0)
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
                     r2 = xi * xi + yi * yi;
                     if ((rminsq <= r2) && (r2 <= rmaxsq)) {
                        if ((yi * cm - xi * sm) < 0)
                           return s;
                     }
                  }
               }
            }
         }
         // We may also cross rmin, (+) solution
         if (rdotn >= 0)
            return TGeoShape::Big();
         if (cdfi >= 0)
            return TGeoShape::Big();
         TGeoTube::DistToTube(rsq, nsq, rdotn, rmin, b, d);
         if (d > 0) {
            s = -b + d;
            if (s > 0) {
               zi = point[2] + s * dir[2];
               if (TMath::Abs(zi) <= dz) {
                  xi = point[0] + s * dir[0];
                  yi = point[1] + s * dir[1];
                  if ((xi * cm + yi * sm) >= rmin * cdfi)
                     return s;
               }
            }
         }
         return TGeoShape::Big();
      }
      // we are on rmin boundary: we may cross again rmin or a phi facette
      if (rdotn >= 0)
         return 0.0;
      TGeoTube::DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d > 0) {
         s = -b + d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               // now check phi range
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               if ((xi * cm + yi * sm) >= rmin * cdfi)
                  return s;
               // now we really have to check any phi crossing
               Double_t un = -dir[0] * s1 + dir[1] * c1;
               if (un > 0) {
                  s = point[0] * s1 - point[1] * c1;
                  if (s >= 0) {
                     s /= un;
                     zi = point[2] + s * dir[2];
                     if (TMath::Abs(zi) <= dz) {
                        xi = point[0] + s * dir[0];
                        yi = point[1] + s * dir[1];
                        r2 = xi * xi + yi * yi;
                        if ((rminsq <= r2) && (r2 <= rmaxsq)) {
                           if ((yi * cm - xi * sm) <= 0) {
                              if (s < snxt)
                                 snxt = s;
                           }
                        }
                     }
                  }
               }
               un = dir[0] * s2 - dir[1] * c2;
               if (un > 0) {
                  s = (point[1] * c2 - point[0] * s2) / un;
                  if (s >= 0 && s < snxt) {
                     zi = point[2] + s * dir[2];
                     if (TMath::Abs(zi) <= dz) {
                        xi = point[0] + s * dir[0];
                        yi = point[1] + s * dir[1];
                        r2 = xi * xi + yi * yi;
                        if ((rminsq <= r2) && (r2 <= rmaxsq)) {
                           if ((yi * cm - xi * sm) >= 0) {
                              return s;
                           }
                        }
                     }
                  }
               }
               return snxt;
            }
         }
      }
      return TGeoShape::Big();
   }
   // only r>rmax has to be considered
   if (TMath::Abs(nsq) < TGeoShape::Tolerance())
      return TGeoShape::Big();
   if (rsq >= rmax * rmax) {
      if (rdotn >= 0)
         return TGeoShape::Big();
      TGeoTube::DistToTube(rsq, nsq, rdotn, rmax, b, d);
      if (d > 0) {
         s = -b - d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               cpsi = xi * cm + yi * sm;
               if (cpsi >= rmax * cdfi)
                  return s;
            }
         }
      }
   }
   // check inner cylinder
   if (rmin > 0) {
      TGeoTube::DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d > 0) {
         s = -b + d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz) {
               xi = point[0] + s * dir[0];
               yi = point[1] + s * dir[1];
               cpsi = xi * cm + yi * sm;
               if (cpsi >= rmin * cdfi)
                  snxt = s;
            }
         }
      }
   }
   // check phi planes
   Double_t un = -dir[0] * s1 + dir[1] * c1;
   if (un > 0) {
      s = point[0] * s1 - point[1] * c1;
      if (s >= 0) {
         s /= un;
         zi = point[2] + s * dir[2];
         if (TMath::Abs(zi) <= dz) {
            xi = point[0] + s * dir[0];
            yi = point[1] + s * dir[1];
            r2 = xi * xi + yi * yi;
            if ((rminsq <= r2) && (r2 <= rmaxsq)) {
               if ((yi * cm - xi * sm) <= 0) {
                  if (s < snxt)
                     snxt = s;
               }
            }
         }
      }
   }
   un = dir[0] * s2 - dir[1] * c2;
   if (un > 0) {
      s = point[1] * c2 - point[0] * s2;
      if (s >= 0) {
         s /= un;
         zi = point[2] + s * dir[2];
         if (TMath::Abs(zi) <= dz) {
            xi = point[0] + s * dir[0];
            yi = point[1] + s * dir[1];
            r2 = xi * xi + yi * yi;
            if ((rminsq <= r2) && (r2 <= rmaxsq)) {
               if ((yi * cm - xi * sm) >= 0) {
                  if (s < snxt)
                     snxt = s;
               }
            }
         }
      }
   }
   return snxt;
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Divide this tube segment shape belonging to volume "voldiv" into ndiv volumes
// /// called divname, from start position with the given step. Returns pointer
// /// to created division cell volume in case of Z divisions. For radialdivision
// /// creates all volumes with different shapes and returns pointer to volume that
// /// was divided. In case a wrong division axis is supplied, returns pointer to
// /// volume that was divided.

TGeoVolume *
TGeoVGTubeSeg::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
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
   case 1: //---                 R division
      finder = new TGeoPatternCylR(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         shape = new TGeoVGTubeSeg(start + id * step, start + (id + 1) * step, z(), sphi() * TMath::RadToDeg(), (sphi() + Base_t::dphi()) * TMath::RadToDeg());
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "R";
         voldiv->AddNodeOffset(vol, id, 0, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 2: //---                 Phi division
      dphi = Base_t::dphi() * TMath::RadToDeg();
      if (dphi < 0)
         dphi += 360.;
      if (step <= 0) {
         step = dphi / ndiv;
         start = sphi() * TMath::RadToDeg();
         end = (sphi() + Base_t::dphi()) * TMath::RadToDeg();
      }
      finder = new TGeoPatternCylPhi(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGTubeSeg(rmin(), rmax(), z(), -step / 2, step / 2);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      opt = "Phi";
      for (id = 0; id < ndiv; id++) {
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 3: //---                  Z division
      finder = new TGeoPatternZ(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGTubeSeg(rmin(), rmax(), step / 2, sphi() * TMath::RadToDeg(), (sphi() + Base_t::dphi())* TMath::RadToDeg());
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      opt = "Z";
      for (id = 0; id < ndiv; id++) {
         voldiv->AddNodeOffset(vol, id, start + step / 2 + id * step, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   default: Error("Divide", "In shape %s wrong axis type for division", GetName()); return nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGTubeSeg::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
      xlo = sphi() * TMath::RadToDeg();
      xhi = (sphi() + dphi()) * TMath::RadToDeg();
      dx = xhi - xlo;
      return dx ;
   case 3:
      xlo = -z();
      xhi = z();
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGTubeSeg::GetBoundingCylinder(Double_t *param) const
{
   param[0] = rmin();
   param[0] *= param[0];
   param[1] = rmax();
   param[1] *= param[1];
   param[2] = sphi() * TMath::RadToDeg();
   param[3] = (sphi() + dphi()) * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGTubeSeg::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
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
   if (z() < 0) {
      // TODO: use typeinfo 
      if (dynamic_cast<TGeoTube*>(mother)) {
        dz = ((TGeoTube *)mother)->GetDz();
      }
      // if (dynamic_cast<TGeoVGTube*>(mother)) {
      //   dz = ((TGeoVGTube *)mother)->GetDz();
      // }
   }
   if (Base_t::rmin() < 0) {
      if (dynamic_cast<TGeoTube*>(mother)) {
        rmin = ((TGeoTube *)mother)->GetRmin();
      }
      // if (dynamic_cast<TGeoVGTube*>(mother)) {
      //   rmin = ((TGeoVGTube *)mother)->GetRmin();
      // }
   }
   if ((Base_t::rmax() < 0) || (Base_t::rmax() <= Base_t::rmin())) {
      if (dynamic_cast<TGeoTube*>(mother)) {
        rmax = ((TGeoTube *)mother)->GetRmax();
      }
      // if (dynamic_cast<TGeoVGTube*>(mother)) {
      //   rmax = ((TGeoVGTube *)mother)->GetRmax();
      // }
   }

   return (new TGeoVGTubeSeg(GetName(), rmin, rmax, dz, sphi() * TMath::RadToDeg(), (sphi() + dphi()) * TMath::RadToDeg()));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGTubeSeg::InspectShape() const
{
   printf("*** Shape %s: TGeoVGTubeSeg ***\n", GetName());
   printf("    Rmin = %11.5f\n", rmin());
   printf("    Rmax = %11.5f\n", rmax());
   printf("    dz   = %11.5f\n", z());
   printf("    phi1 = %11.5f\n", sphi()* TMath::RadToDeg());
   printf("    phi2 = %11.5f\n", (sphi() + dphi()) * TMath::RadToDeg());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGTubeSeg::MakeBuffer3D() const
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

void TGeoVGTubeSeg::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t c = GetBasicColor();

   memset(buff.fSegs, 0, buff.NbSegs() * 3 * sizeof(Int_t));
   for (i = 0; i < 4; i++) {
      for (j = 1; j < n; j++) {
         buff.fSegs[(i * n + j - 1) * 3] = c;
         buff.fSegs[(i * n + j - 1) * 3 + 1] = i * n + j - 1;
         buff.fSegs[(i * n + j - 1) * 3 + 2] = i * n + j;
      }
   }
   for (i = 4; i < 6; i++) {
      for (j = 0; j < n; j++) {
         buff.fSegs[(i * n + j) * 3] = c + 1;
         buff.fSegs[(i * n + j) * 3 + 1] = (i - 4) * n + j;
         buff.fSegs[(i * n + j) * 3 + 2] = (i - 2) * n + j;
      }
   }
   for (i = 6; i < 8; i++) {
      for (j = 0; j < n; j++) {
         buff.fSegs[(i * n + j) * 3] = c;
         buff.fSegs[(i * n + j) * 3 + 1] = 2 * (i - 6) * n + j;
         buff.fSegs[(i * n + j) * 3 + 2] = (2 * (i - 6) + 1) * n + j;
      }
   }

   Int_t indx = 0;
   memset(buff.fPols, 0, buff.NbPols() * 6 * sizeof(Int_t));
   i = 0;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (4 + i) * n + j + 1;
      buff.fPols[indx++] = (2 + i) * n + j;
      buff.fPols[indx++] = (4 + i) * n + j;
      buff.fPols[indx++] = i * n + j;
   }
   i = 1;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = i * n + j;
      buff.fPols[indx++] = (4 + i) * n + j;
      buff.fPols[indx++] = (2 + i) * n + j;
      buff.fPols[indx++] = (4 + i) * n + j + 1;
   }
   i = 2;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c + i;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (i - 2) * 2 * n + j;
      buff.fPols[indx++] = (4 + i) * n + j;
      buff.fPols[indx++] = ((i - 2) * 2 + 1) * n + j;
      buff.fPols[indx++] = (4 + i) * n + j + 1;
   }
   i = 3;
   for (j = 0; j < n - 1; j++) {
      buff.fPols[indx++] = c + i;
      buff.fPols[indx++] = 4;
      buff.fPols[indx++] = (4 + i) * n + j + 1;
      buff.fPols[indx++] = ((i - 2) * 2 + 1) * n + j;
      buff.fPols[indx++] = (4 + i) * n + j;
      buff.fPols[indx++] = (i - 2) * 2 * n + j;
   }
   buff.fPols[indx++] = c + 2;
   buff.fPols[indx++] = 4;
   buff.fPols[indx++] = 6 * n;
   buff.fPols[indx++] = 4 * n;
   buff.fPols[indx++] = 7 * n;
   buff.fPols[indx++] = 5 * n;
   buff.fPols[indx++] = c + 2;
   buff.fPols[indx++] = 4;
   buff.fPols[indx++] = 6 * n - 1;
   buff.fPols[indx++] = 8 * n - 1;
   buff.fPols[indx++] = 5 * n - 1;
   buff.fPols[indx++] = 7 * n - 1;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: method to compute the closest distance from given point to this shape.

Double_t TGeoVGTubeSeg::SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz,
                              Double_t phi1d, Double_t phi2d, Int_t skipz)
{
   Double_t saf[3];
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t r = TMath::Sqrt(rsq);

   switch (skipz) {
   case 1: // skip lower Z plane
      saf[0] = dz - point[2];
      break;
   case 2: // skip upper Z plane
      saf[0] = dz + point[2];
      break;
   case 3: // skip both
      saf[0] = TGeoShape::Big();
      break;
   default: saf[0] = dz - TMath::Abs(point[2]);
   }

   if (in) {
      saf[1] = r - rmin;
      saf[2] = rmax - r;
      Double_t safe = saf[TMath::LocMin(3, saf)];
      if ((phi2d - phi1d) >= 360.)
         return safe;
      Double_t safphi = TGeoShape::SafetyPhi(point, in, phi1d, phi2d);
      return TMath::Min(safe, safphi);
   }
   // Point expected to be outside
   saf[0] = -saf[0];
   Bool_t inphi = kFALSE;
   Double_t phi1 = phi1d * TMath::DegToRad();
   Double_t phi2 = phi2d * TMath::DegToRad();

   Double_t fio = 0.5 * (phi1 + phi2);
   Double_t cm = TMath::Cos(fio);
   Double_t sm = TMath::Sin(fio);
   Double_t cpsi = point[0] * cm + point[1] * sm;
   Double_t dfi = 0.5 * (phi2 - phi1);
   Double_t cdfi = TMath::Cos(dfi);
   if (cpsi > r * cdfi - TGeoShape::Tolerance())
      inphi = kTRUE;
   if (inphi) {
      saf[1] = rmin - r;
      saf[2] = r - rmax;
      Double_t safe = saf[TMath::LocMax(3, saf)];
      safe = TMath::Max(0., safe);
      return safe;
   }
   // Point outside the phi range
   // Compute projected radius of the (r,phi) position vector onto
   // phi1 and phi2 edges and take the maximum for choosing the side.
   Double_t c1 = TMath::Cos(phi1);
   Double_t s1 = TMath::Sin(phi1);
   Double_t c2 = TMath::Cos(phi2);
   Double_t s2 = TMath::Sin(phi2);

   Double_t rproj = TMath::Max(point[0] * c1 + point[1] * s1, point[0] * c2 + point[1] * s2);
   saf[1] = rmin - rproj;
   saf[2] = rproj - rmax;
   Double_t safe = TMath::Max(saf[1], saf[2]);
   if ((phi2d - phi1d) >= 360.)
      return TMath::Max(safe, saf[0]);
   if (safe > 0) {
      // rproj not within (rmin,rmax) - > no need to calculate safphi
      safe = TMath::Sqrt(rsq - rproj * rproj + safe * safe);
      return (saf[0] < 0) ? safe : TMath::Sqrt(safe * safe + saf[0] * saf[0]);
   }
   Double_t safphi = TGeoShape::SafetyPhi(point, in, phi1d, phi2d);
   return (saf[0] < 0) ? safphi : TMath::Sqrt(saf[0] * saf[0] + safphi * safphi);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGTubeSeg::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << rmin() << ";" << std::endl;
   out << "   rmax = " << rmax() << ";" << std::endl;
   out << "   dz   = " << z() << ";" << std::endl;
   out << "   phi1 = " << sphi() * TMath::RadToDeg() << ";" << std::endl;
   out << "   phi2 = " << (sphi() + dphi()) * TMath::RadToDeg() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGTubeSeg(\"" << GetName() << "\",rmin,rmax,dz,phi1,phi2);"
       << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi1

Double_t TGeoVGTubeSeg::GetPhi1() const
{
   return sphi() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi2

Double_t TGeoVGTubeSeg::GetPhi2() const
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
/// Set dimensions of the tube segment starting from a list.

void TGeoVGTubeSeg::SetDimensions(Double_t *param)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t dz = param[2];
   Double_t phi1 = param[3];
   Double_t phi2 = param[4];
   SetTubsDimensions(rmin, rmax, dz, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoVGTubeSeg::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints / 2) * 2) {
      Error("GetPointsOnSegments", "Npoints must be even number");
      return kFALSE;
   }
   Int_t nc = (Int_t)TMath::Sqrt(0.5 * npoints);
   Double_t dphi = (Base_t::dphi()) / (nc - 1);
   Double_t phi = 0;
   Double_t phi1 = sphi();
   Int_t ntop = npoints / 2 - nc * (nc - 1);
   Double_t dz = 2 * z() / (nc - 1);
   Double_t z = 0;
   Int_t icrt = 0;
   Int_t nphi = nc;
   // loop z sections
   for (Int_t i = 0; i < nc; i++) {
      if (i == (nc - 1)) {
         nphi = ntop;
         dphi = (Base_t::dphi()) / (nphi - 1);
      }
      z = -Base_t::z() + i * dz;
      // loop points on circle sections
      for (Int_t j = 0; j < nphi; j++) {
         phi = phi1 + j * dphi;
         array[icrt++] = rmin() * TMath::Cos(phi);
         array[icrt++] = rmin() * TMath::Sin(phi);
         array[icrt++] = z;
         array[icrt++] = rmax() * TMath::Cos(phi);
         array[icrt++] = rmax() * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Create tube segment mesh points.

void TGeoVGTubeSeg::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = sphi() ;
   phi2 = sphi() + Base_t::dphi();
   if (phi2 < phi1)
      phi2 += TMath::TwoPi();
   n = gGeoManager->GetNsegments() + 1;

   dphi = (phi2 - phi1) / (n - 1);
   dz = z();

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi);
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi);
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create tube segment mesh points.

void TGeoVGTubeSeg::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   Double_t phi, phi1, phi2, dphi;
   phi1 = sphi();
   phi2 = sphi() + Base_t::dphi();
   if (phi2 < phi1)
      phi2 += TMath::TwoPi();
   n = gGeoManager->GetNsegments() + 1;

   dphi = (phi2 - phi1) / (n - 1);
   dz = z();

   if (points) {
      Int_t indx = 0;

      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi);
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = (phi1 + j * dphi);
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
         indx++;
         points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGTubeSeg::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   nvert = n * 4;
   nsegs = n * 8;
   npols = n * 4 - 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGTubeSeg::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t numPoints = n * 4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGTubeSeg::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGTubeSeg::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DTubeSeg buffer;
   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      // These from TBuffer3DTube / TGeoTube
      buffer.fRadiusInner = rmin();
      buffer.fRadiusOuter = rmax();
      buffer.fHalfLength = z();
      buffer.fPhiMin = sphi() * TMath::RadToDeg();
      buffer.fPhiMax = (sphi() + dphi()) * TMath::RadToDeg();
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

template class TGeoVGAdapter<vecgeom::cxx::SUnplacedTube<vecgeom::cxx::TubeTypes::UniversalTube>>;
