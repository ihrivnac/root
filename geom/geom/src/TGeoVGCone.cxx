// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoVGCone::Contains() and DistFromInside() implemented by Mihaela Gheata

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

// ClassImp(TGeoVGCone);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGCone::TGeoVGCone()
 : Base_t("", 0., 0., 0., 0., 0., 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoCone);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGCone::TGeoVGCone(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2)
 : Base_t("", rmin1, rmax1, rmin2, rmax2, dz, 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoCone);
   // SetConeDimensions(dz, rmin1, rmax1, rmin2, rmax2);
   if ((dz < 0) || (rmin1 < 0) || (rmax1 < 0) || (rmin2 < 0) || (rmax2 < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGCone::TGeoVGCone(const char *name, Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2)
 : Base_t(name, rmin1, rmax1, rmin2, rmax2, dz, 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoCone);
   // SetConeDimensions(dz, rmin1, rmax1, rmin2, rmax2);
   if ((dz < 0) || (rmin1 < 0) || (rmax1 < 0) || (rmin2 < 0) || (rmax2 < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = dz
///  - param[1] = Rmin1
///  - param[2] = Rmax1
///  - param[3] = Rmin2
///  - param[4] = Rmax2

TGeoVGCone::TGeoVGCone(Double_t *param)
 : Base_t("", param[1], param[2], param[3], param[3], param[0], 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoCone);
   // SetDimensions(param);
   if ((GetDz() < 0) || (GetRmin1() < 0) || (GetRmax1() < 0) || (GetRmin2() < 0) || (GetRmax2() < 0))
      SetShapeBit(kGeoRunTimeShape);
   else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Computes capacity of the shape in [length^3]

Double_t TGeoVGCone::Capacity(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2)
{
   Double_t capacity = (2. * dz * TMath::Pi() / 3.) *
                       (rmax1 * rmax1 + rmax2 * rmax2 + rmax1 * rmax2 - rmin1 * rmin1 - rmin2 * rmin2 - rmin1 * rmin2);
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGCone::~TGeoVGCone() {}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the sphere

void TGeoVGCone::ComputeBBox()
{
   TGeoBBox *box = (TGeoBBox *)this;
   box->SetBoxDimensions(TMath::Max(GetRmax1(), GetRmax2()), TMath::Max(GetRmax1(), GetRmax2()), GetDz());
   memset(fOrigin, 0, 3 * sizeof(Double_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute normal to closest surface from POINT.

void TGeoVGCone::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t dz, Double_t rmin1,
                              Double_t rmax1, Double_t rmin2, Double_t rmax2)
{
   Double_t safe, phi;
   memset(norm, 0, 3 * sizeof(Double_t));
   phi = TMath::ATan2(point[1], point[0]);
   Double_t cphi = TMath::Cos(phi);
   Double_t sphi = TMath::Sin(phi);
   Double_t ro1 = 0.5 * (rmin1 + rmin2);
   Double_t tg1 = 0.5 * (rmin2 - rmin1) / dz;
   Double_t cr1 = 1. / TMath::Sqrt(1. + tg1 * tg1);
   Double_t ro2 = 0.5 * (rmax1 + rmax2);
   Double_t tg2 = 0.5 * (rmax2 - rmax1) / dz;
   Double_t cr2 = 1. / TMath::Sqrt(1. + tg2 * tg2);

   Double_t r = TMath::Sqrt(point[0] * point[0] + point[1] * point[1]);
   Double_t rin = tg1 * point[2] + ro1;
   Double_t rout = tg2 * point[2] + ro2;
   safe = (ro1 > 0) ? (TMath::Abs((r - rin) * cr1)) : TGeoShape::Big();
   norm[0] = cr1 * cphi;
   norm[1] = cr1 * sphi;
   norm[2] = -tg1 * cr1;
   if (TMath::Abs((rout - r) * cr2) < safe) {
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
/// Static: Compute distance from inside point to surface of the cone (static)
/// Boundary safe algorithm.

Double_t TGeoVGCone::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                   Double_t rmax1, Double_t rmin2, Double_t rmax2)
{
   if (dz <= 0)
      return TGeoShape::Big();
   // compute distance to surface
   // Do Z
   Double_t sz = TGeoShape::Big();
   if (dir[2]) {
      sz = (TMath::Sign(dz, dir[2]) - point[2]) / dir[2];
      if (sz <= 0)
         return 0.0;
   }
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t zinv = 1. / dz;
   Double_t rin = 0.5 * (rmin1 + rmin2 + (rmin2 - rmin1) * point[2] * zinv);
   // Do Rmin
   Double_t b, delta, zi;
   if (rin > 0) {
      // Protection in case point is actually outside the cone
      if (rsq < rin * (rin + TGeoShape::Tolerance())) {
         Double_t ddotn =
            point[0] * dir[0] + point[1] * dir[1] + 0.5 * (rmin1 - rmin2) * dir[2] * zinv * TMath::Sqrt(rsq);
         if (ddotn <= 0)
            return 0.0;
      } else {
         TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
         if (delta > 0) {
            Double_t sr = -b - delta;
            if (sr > 0) {
               zi = point[2] + sr * dir[2];
               if (TMath::Abs(zi) <= dz)
                  return TMath::Min(sz, sr);
            }
            sr = -b + delta;
            if (sr > 0) {
               zi = point[2] + sr * dir[2];
               if (TMath::Abs(zi) <= dz)
                  return TMath::Min(sz, sr);
            }
         }
      }
   }
   // Do Rmax
   Double_t rout = 0.5 * (rmax1 + rmax2 + (rmax2 - rmax1) * point[2] * zinv);
   if (rsq > rout * (rout - TGeoShape::Tolerance())) {
      Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] + 0.5 * (rmax1 - rmax2) * dir[2] * zinv * TMath::Sqrt(rsq);
      if (ddotn >= 0)
         return 0.0;
      TGeoVGCone::DistToCone(point, dir, dz, rmax1, rmax2, b, delta);
      if (delta < 0)
         return 0.0;
      Double_t sr = -b + delta;
      if (sr < 0)
         return sz;
      if (TMath::Abs(-b - delta) > sr)
         return sz;
      zi = point[2] + sr * dir[2];
      if (TMath::Abs(zi) <= dz)
         return TMath::Min(sz, sr);
      return sz;
   }
   TGeoVGCone::DistToCone(point, dir, dz, rmax1, rmax2, b, delta);
   if (delta > 0) {
      Double_t sr = -b - delta;
      if (sr > 0) {
         zi = point[2] + sr * dir[2];
         if (TMath::Abs(zi) <= dz)
            return TMath::Min(sz, sr);
      }
      sr = -b + delta;
      if (sr > TGeoShape::Tolerance()) {
         zi = point[2] + sr * dir[2];
         if (TMath::Abs(zi) <= dz)
            return TMath::Min(sz, sr);
      }
   }
   return sz;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute distance from outside point to surface of the tube
/// Boundary safe algorithm.

Double_t TGeoVGCone::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t dz, Double_t rmin1,
                                    Double_t rmax1, Double_t rmin2, Double_t rmax2)
{
   // compute distance to Z planes
   if (dz <= 0)
      return TGeoShape::Big();
   Double_t snxt;
   Double_t xp, yp, zp;
   Bool_t inz = kTRUE;

   if (point[2] <= -dz) {
      if (dir[2] <= 0)
         return TGeoShape::Big();
      snxt = (-dz - point[2]) / dir[2];
      xp = point[0] + snxt * dir[0];
      yp = point[1] + snxt * dir[1];
      Double_t r2 = xp * xp + yp * yp;
      if ((r2 >= rmin1 * rmin1) && (r2 <= rmax1 * rmax1))
         return snxt;
      inz = kFALSE;
   } else {
      if (point[2] >= dz) {
         if (dir[2] >= 0)
            return TGeoShape::Big();
         snxt = (dz - point[2]) / dir[2];
         xp = point[0] + snxt * dir[0];
         yp = point[1] + snxt * dir[1];
         Double_t r2 = xp * xp + yp * yp;
         if ((r2 >= rmin2 * rmin2) && (r2 <= rmax2 * rmax2))
            return snxt;
         inz = kFALSE;
      }
   }

   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t dzinv = 1. / dz;
   Double_t ro1 = 0.5 * (rmin1 + rmin2);
   Bool_t hasrmin = (ro1 > 0) ? kTRUE : kFALSE;
   Double_t tg1 = 0.;
   Double_t rin = 0.;
   Bool_t inrmin = kTRUE; // r>=rmin
   if (hasrmin) {
      tg1 = 0.5 * (rmin2 - rmin1) * dzinv;
      rin = ro1 + tg1 * point[2];
      if (rin > 0 && rsq < rin * (rin - TGeoShape::Tolerance()))
         inrmin = kFALSE;
   }
   Double_t ro2 = 0.5 * (rmax1 + rmax2);
   Double_t tg2 = 0.5 * (rmax2 - rmax1) * dzinv;
   Double_t rout = tg2 * point[2] + ro2;
   Bool_t inrmax = kFALSE;
   if (rout > 0 && rsq < rout * (rout + TGeoShape::Tolerance()))
      inrmax = kTRUE;
   Bool_t in = inz & inrmin & inrmax;
   Double_t b, delta;
   // If inside cone, we are most likely on a boundary within machine precision.
   if (in) {
      Double_t r = TMath::Sqrt(rsq);
      Double_t safz = dz - TMath::Abs(point[2]); // positive
      Double_t safrmin = (hasrmin) ? (r - rin) : TGeoShape::Big();
      Double_t safrmax = rout - r;
      if (safz <= safrmin && safz <= safrmax) {
         // on Z boundary
         if (point[2] * dir[2] < 0)
            return 0.0;
         return TGeoShape::Big();
      }
      if (safrmax < safrmin) {
         // on rmax boundary
         Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] - tg2 * dir[2] * r;
         if (ddotn <= 0)
            return 0.0;
         return TGeoShape::Big();
      }
      // on rmin boundary
      Double_t ddotn = point[0] * dir[0] + point[1] * dir[1] - tg1 * dir[2] * r;
      if (ddotn >= 0)
         return 0.0;
      // we can cross (+) solution of rmin
      TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);

      if (delta < 0)
         return 0.0;
      snxt = -b + delta;
      if (snxt < 0)
         return TGeoShape::Big();
      if (TMath::Abs(-b - delta) > snxt)
         return TGeoShape::Big();
      zp = point[2] + snxt * dir[2];
      if (TMath::Abs(zp) <= dz)
         return snxt;
      return TGeoShape::Big();
   }

   // compute distance to inner cone
   snxt = TGeoShape::Big();
   if (!inrmin) {
      // ray can cross inner cone (but not only!)
      TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
      if (delta < 0)
         return TGeoShape::Big();
      snxt = -b + delta;
      if (snxt > 0) {
         zp = point[2] + snxt * dir[2];
         if (TMath::Abs(zp) <= dz)
            return snxt;
      }
      snxt = -b - delta;
      if (snxt > 0) {
         zp = point[2] + snxt * dir[2];
         if (TMath::Abs(zp) <= dz)
            return snxt;
      }
      snxt = TGeoShape::Big();
   } else {
      if (hasrmin) {
         TGeoVGCone::DistToCone(point, dir, dz, rmin1, rmin2, b, delta);
         if (delta > 0) {
            Double_t din = -b + delta;
            if (din > 0) {
               zp = point[2] + din * dir[2];
               if (TMath::Abs(zp) <= dz)
                  snxt = din;
            }
         }
      }
   }

   if (inrmax)
      return snxt;
   // We can cross outer cone, both solutions possible
   // compute distance to outer cone
   TGeoVGCone::DistToCone(point, dir, dz, rmax1, rmax2, b, delta);
   if (delta < 0)
      return snxt;
   Double_t dout = -b - delta;
   if (dout > 0 && dout < snxt) {
      zp = point[2] + dout * dir[2];
      if (TMath::Abs(zp) <= dz)
         return dout;
   }
   dout = -b + delta;
   if (dout <= 0 || dout > snxt)
      return snxt;
   zp = point[2] + dout * dir[2];
   if (TMath::Abs(zp) <= dz)
      return dout;
   return snxt;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute distance to a conical surface with :
/// - r1, z1 - radius and Z position of lower base
/// - r2, z2 - radius and Z position of upper base

void TGeoVGCone::DistToCone(const Double_t *point, const Double_t *dir, Double_t dz, Double_t r1, Double_t r2,
                          Double_t &b, Double_t &delta)
{
   b = 0;
   delta = -1.;
   if (dz < 0)
      return;
   Double_t ro0 = 0.5 * (r1 + r2);
   Double_t tz = 0.5 * (r2 - r1) / dz;
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t rc = ro0 + point[2] * tz;

   Double_t a = dir[0] * dir[0] + dir[1] * dir[1] - tz * tz * dir[2] * dir[2];
   b = point[0] * dir[0] + point[1] * dir[1] - tz * rc * dir[2];
   Double_t c = rsq - rc * rc;

   if (TMath::Abs(a) < TGeoShape::Tolerance()) {
      if (TMath::Abs(b) < TGeoShape::Tolerance())
         return;
      b = 0.5 * c / b;
      delta = 0.;
      return;
   }
   a = 1. / a;
   b *= a;
   c *= a;
   delta = b * b - c;
   if (delta > 0) {
      delta = TMath::Sqrt(delta);
   } else {
      delta = -1.;
   }
}

// ////////////////////////////////////////////////////////////////////////////////
// /// compute closest distance from point px,py to each corner

// Int_t TGeoVGCone::DistancetoPrimitive(Int_t px, Int_t py)
// {
//    Int_t n = gGeoManager->GetNsegments();
//    const Int_t numPoints = 4 * n;
//    return ShapeDistancetoPrimitive(numPoints, px, py);
// }

////////////////////////////////////////////////////////////////////////////////
/// Divide this cone shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. For Z division
/// creates all volumes with different shapes and returns pointer to volume that
/// was divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *
TGeoVGCone::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Int_t id;
   Double_t end = start + ndiv * step;
   switch (iaxis) {
   case 1: //---              R division
      Error("Divide", "division of a cone on R not implemented");
      return nullptr;
   case 2: // ---             Phi division
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
   case 3: //---               Z division
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      finder = new TGeoPatternZ(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         Double_t z1 = start + id * step;
         Double_t z2 = start + (id + 1) * step;
         Double_t rmin1n = 0.5 * (GetRmin1() * (GetDz() - z1) + GetRmin2() * (GetDz() + z1)) / GetDz();
         Double_t rmax1n = 0.5 * (GetRmax1() * (GetDz() - z1) + GetRmax2() * (GetDz() + z1)) / GetDz();
         Double_t rmin2n = 0.5 * (GetRmin1() * (GetDz() - z2) + GetRmin2() * (GetDz() + z2)) / GetDz();
         Double_t rmax2n = 0.5 * (GetRmax1() * (GetDz() - z2) + GetRmax2() * (GetDz() + z2)) / GetDz();
         shape = new TGeoVGCone(0.5 * step, rmin1n, rmax1n, rmin2n, rmax2n);
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
/// Returns name of axis IAXIS.

const char *TGeoVGCone::GetAxisName(Int_t iaxis) const
{
   switch (iaxis) {
   case 1: return "R";
   case 2: return "PHI";
   case 3: return "Z";
   default: return "undefined";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGCone::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 2:
      xlo = 0.;
      xhi = 360.;
      return 360.;
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
/// is the following : Rmin, Rmax, Phi1, Phi2, dZ

void TGeoVGCone::GetBoundingCylinder(Double_t *param) const
{
   param[0] = TMath::Min(GetRmin1(), GetRmin2()); // Rmin
   param[0] *= param[0];
   param[1] = TMath::Max(GetRmax1(), GetRmax2()); // Rmax
   param[1] *= param[1];
   param[2] = 0.;   // Phi1
   param[3] = 360.; // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGCone::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (!mother->TestShapeBit(kGeoCone)) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t rmin1, rmax1, rmin2, rmax2, dz;
   rmin1 = GetRmin1();
   rmax1 = GetRmax1();
   rmin2 = GetRmin2();
   rmax2 = GetRmax2();
   dz = Base_t::GetDz();
   if (Base_t::GetDz() < 0)
      dz = ((TGeoVGCone *)mother)->GetDz();
   if (GetRmin1() < 0)
      rmin1 = ((TGeoVGCone *)mother)->GetRmin1();
   if (GetRmax1() < 0)
      rmax1 = ((TGeoVGCone *)mother)->GetRmax1();
   if (GetRmin2() < 0)
      rmin2 = ((TGeoVGCone *)mother)->GetRmin2();
   if (GetRmax2() < 0)
      rmax2 = ((TGeoVGCone *)mother)->GetRmax2();

   return (new TGeoVGCone(GetName(), dz, rmin1, rmax1, rmin2, rmax2));
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoVGCone::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints / 2) * 2) {
      Error("GetPointsOnSegments", "Npoints must be even number");
      return kFALSE;
   }
   Bool_t hasrmin = (GetRmin1() > 0 || GetRmin2() > 0) ? kTRUE : kFALSE;
   Int_t nc = 0;
   if (hasrmin)
      nc = (Int_t)TMath::Sqrt(0.5 * npoints);
   else
      nc = (Int_t)TMath::Sqrt(1. * npoints);
   Double_t dphi = TMath::TwoPi() / nc;
   Double_t phi = 0;
   Int_t ntop = 0;
   if (hasrmin)
      ntop = npoints / 2 - nc * (nc - 1);
   else
      ntop = npoints - nc * (nc - 1);
   Double_t dz = 2 * Base_t::GetDz() / (nc - 1);
   Double_t z = 0;
   Int_t icrt = 0;
   Int_t nphi = nc;
   Double_t rmin = 0.;
   Double_t rmax = 0.;
   // loop z sections
   for (Int_t i = 0; i < nc; i++) {
      if (i == (nc - 1))
         nphi = ntop;
      z = -Base_t::GetDz() + i * dz;
      if (hasrmin)
         rmin = 0.5 * (GetRmin1() + GetRmin2()) + 0.5 * (GetRmin2() - GetRmin1()) * z / Base_t::GetDz();
      rmax = 0.5 * (GetRmax1() + GetRmax2()) + 0.5 * (GetRmax2() - GetRmax1()) * z / Base_t::GetDz();
      // loop points on circle sections
      for (Int_t j = 0; j < nphi; j++) {
         phi = j * dphi;
         if (hasrmin) {
            array[icrt++] = rmin * TMath::Cos(phi);
            array[icrt++] = rmin * TMath::Sin(phi);
            array[icrt++] = z;
         }
         array[icrt++] = rmax * TMath::Cos(phi);
         array[icrt++] = rmax * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGCone::InspectShape() const
{
   printf("*** Shape %s TGeoVGCone ***\n", GetName());
   printf("    dz    =: %11.5f\n", GetDz());
   printf("    Rmin1 = %11.5f\n", GetRmin1());
   printf("    Rmax1 = %11.5f\n", GetRmax1());
   printf("    Rmin2 = %11.5f\n", GetRmin2());
   printf("    Rmax2 = %11.5f\n", GetRmax2());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGCone::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t nbPnts = 4 * n;
   Int_t nbSegs = 8 * n;
   Int_t nbPols = 4 * n;
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

void TGeoVGCone::SetSegsAndPols(TBuffer3D &buffer) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments();
   Int_t c = GetBasicColor();

   for (i = 0; i < 4; i++) {
      for (j = 0; j < n; j++) {
         buffer.fSegs[(i * n + j) * 3] = c;
         buffer.fSegs[(i * n + j) * 3 + 1] = i * n + j;
         buffer.fSegs[(i * n + j) * 3 + 2] = i * n + j + 1;
      }
      buffer.fSegs[(i * n + j - 1) * 3 + 2] = i * n;
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
   i = 0;
   for (j = 0; j < n; j++) {
      indx = 6 * (i * n + j);
      buffer.fPols[indx] = c;
      buffer.fPols[indx + 1] = 4;
      buffer.fPols[indx + 5] = i * n + j;
      buffer.fPols[indx + 4] = (4 + i) * n + j;
      buffer.fPols[indx + 3] = (2 + i) * n + j;
      buffer.fPols[indx + 2] = (4 + i) * n + j + 1;
   }
   buffer.fPols[indx + 2] = (4 + i) * n;
   i = 1;
   for (j = 0; j < n; j++) {
      indx = 6 * (i * n + j);
      buffer.fPols[indx] = c;
      buffer.fPols[indx + 1] = 4;
      buffer.fPols[indx + 2] = i * n + j;
      buffer.fPols[indx + 3] = (4 + i) * n + j;
      buffer.fPols[indx + 4] = (2 + i) * n + j;
      buffer.fPols[indx + 5] = (4 + i) * n + j + 1;
   }
   buffer.fPols[indx + 5] = (4 + i) * n;
   i = 2;
   for (j = 0; j < n; j++) {
      indx = 6 * (i * n + j);
      buffer.fPols[indx] = c + i;
      buffer.fPols[indx + 1] = 4;
      buffer.fPols[indx + 2] = (i - 2) * 2 * n + j;
      buffer.fPols[indx + 3] = (4 + i) * n + j;
      buffer.fPols[indx + 4] = ((i - 2) * 2 + 1) * n + j;
      buffer.fPols[indx + 5] = (4 + i) * n + j + 1;
   }
   buffer.fPols[indx + 5] = (4 + i) * n;
   i = 3;
   for (j = 0; j < n; j++) {
      indx = 6 * (i * n + j);
      buffer.fPols[indx] = c + i;
      buffer.fPols[indx + 1] = 4;
      buffer.fPols[indx + 5] = (i - 2) * 2 * n + j;
      buffer.fPols[indx + 4] = (4 + i) * n + j;
      buffer.fPols[indx + 3] = ((i - 2) * 2 + 1) * n + j;
      buffer.fPols[indx + 2] = (4 + i) * n + j + 1;
   }
   buffer.fPols[indx + 2] = (4 + i) * n;
}

////////////////////////////////////////////////////////////////////////////////
/// computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoVGCone::SafetyS(const Double_t *point, Bool_t in, Double_t dz, Double_t rmin1, Double_t rmax1,
                           Double_t rmin2, Double_t rmax2, Int_t skipz)
{
   Double_t saf[4];
   Double_t r = TMath::Sqrt(point[0] * point[0] + point[1] * point[1]);
   //   Double_t rin = tg1*point[2]+ro1;
   //   Double_t rout = tg2*point[2]+ro2;
   switch (skipz) {
   case 1: // skip lower Z plane
      saf[0] = TGeoShape::Big();
      saf[1] = TGeoShape::SafetySeg(r, point[2], rmax2, dz, rmin2, dz, !in);
      break;
   case 2: // skip upper Z plane
      saf[0] = TGeoShape::SafetySeg(r, point[2], rmin1, -dz, rmax1, -dz, !in);
      saf[1] = TGeoShape::Big();
      break;
   case 3: // skip both
      saf[0] = saf[1] = TGeoShape::Big();
      break;
   default:
      saf[0] = TGeoShape::SafetySeg(r, point[2], rmin1, -dz, rmax1, -dz, !in);
      saf[1] = TGeoShape::SafetySeg(r, point[2], rmax2, dz, rmin2, dz, !in);
   }
   // Safety to inner part
   if (rmin1 > 0 || rmin2 > 0)
      saf[2] = TGeoShape::SafetySeg(r, point[2], rmin2, dz, rmin1, -dz, !in);
   else
      saf[2] = TGeoShape::Big();
   saf[3] = TGeoShape::SafetySeg(r, point[2], rmax1, -dz, rmax2, dz, !in);
   return saf[TMath::LocMin(4, saf)];
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGCone::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz    = " << GetDz() << ";" << std::endl;
   out << "   rmin1 = " << GetRmin1() << ";" << std::endl;
   out << "   rmax1 = " << GetRmax1() << ";" << std::endl;
   out << "   rmin2 = " << GetRmin2() << ";" << std::endl;
   out << "   rmax2 = " << GetRmax2() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGCone(\"" << GetName()
       << "\", dz,rmin1,rmax1,rmin2,rmax2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi1

Double_t TGeoVGConeSeg::GetPhi1() const
{
   return GetSPhi() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi2

Double_t TGeoVGConeSeg::GetPhi2() const
{
   return (GetSPhi() + GetDPhi()) * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Set cone dimensions.

void TGeoVGCone::SetConeDimensions(Double_t dz, Double_t rmin1, Double_t rmax1, Double_t rmin2, Double_t rmax2)
{
   if (rmin1 >= 0) {
      if (rmax1 > 0) {
         if (rmin1 <= rmax1) {
            // normal rmin/rmax
            SetRmin1(rmin1);
            SetRmax1(rmax1);
         } else {
            SetRmin1(rmax1);
            SetRmax1(rmin1);
            Warning("SetConeDimensions", "rmin1>rmax1 Switch rmin1<->rmax1");
            SetShapeBit(TGeoShape::kGeoBad);
         }
      } else {
         // run-time
         SetRmin1(rmin1);
         SetRmax1(rmax1);
      }
   } else {
      // run-time
      SetRmin1(rmin1);
      SetRmax1(rmax1);
   }
   if (rmin2 >= 0) {
      if (rmax2 > 0) {
         if (rmin2 <= rmax2) {
            // normal rmin/rmax
            SetRmin2(rmin2);
            SetRmax2(rmax2);
         } else {
            SetRmin2(rmax2);
            SetRmax2(rmin2);
            Warning("SetConeDimensions", "rmin2>rmax2 Switch rmin2<->rmax2");
            SetShapeBit(TGeoShape::kGeoBad);
         }
      } else {
         // run-time
         SetRmin2(rmin2);
         SetRmax2(rmax2);
      }
   } else {
      // run-time
      SetRmin2(rmin2);
      SetRmax2(rmax2);
   }

   SetDz(dz);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Set cone dimensions from an array.

void TGeoVGCone::SetDimensions(Double_t *param)
{
   Double_t dz = param[0];
   Double_t rmin1 = param[1];
   Double_t rmax1 = param[2];
   Double_t rmin2 = param[3];
   Double_t rmax2 = param[4];
   SetConeDimensions(dz, rmin1, rmax1, rmin2, rmax2);
}

////////////////////////////////////////////////////////////////////////////////
/// Create cone mesh points.

void TGeoVGCone::SetPoints(Double_t *points) const
{
   Double_t dz, phi, dphi;
   Int_t j, n;

   n = gGeoManager->GetNsegments();
   dphi = 360. / n;
   dz = Base_t::GetDz();
   Int_t indx = 0;

   if (points) {
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmin1() * TMath::Cos(phi);
         points[indx++] = GetRmin1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmax1() * TMath::Cos(phi);
         points[indx++] = GetRmax1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmin2() * TMath::Cos(phi);
         points[indx++] = GetRmin2() * TMath::Sin(phi);
         points[indx++] = dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmax2() * TMath::Cos(phi);
         points[indx++] = GetRmax2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create cone mesh points.

void TGeoVGCone::SetPoints(Float_t *points) const
{
   Double_t dz, phi, dphi;
   Int_t j, n;

   n = gGeoManager->GetNsegments();
   dphi = 360. / n;
   dz = Base_t::GetDz();
   Int_t indx = 0;

   if (points) {
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmin1() * TMath::Cos(phi);
         points[indx++] = GetRmin1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmax1() * TMath::Cos(phi);
         points[indx++] = GetRmax1() * TMath::Sin(phi);
         points[indx++] = -dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmin2() * TMath::Cos(phi);
         points[indx++] = GetRmin2() * TMath::Sin(phi);
         points[indx++] = dz;
      }

      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = GetRmax2() * TMath::Cos(phi);
         points[indx++] = GetRmax2() * TMath::Sin(phi);
         points[indx++] = dz;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGCone::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments();
   nvert = n * 4;
   nsegs = n * 8;
   npols = n * 4;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGCone::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = n * 4;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill size of this 3-D object

void TGeoVGCone::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGCone::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments();
      Int_t nbPnts = 4 * n;
      Int_t nbSegs = 8 * n;
      Int_t nbPols = 4 * n;
      if (buffer.SetRawSizes(nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }

   // TODO: Can we push this as common down to TGeoShape?
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

#endif
