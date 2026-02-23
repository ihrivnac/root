// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoVGSphere::Contains() DistFromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TGeoSphere.h"
#include "TGeoVGSphere.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>

#include "TGeoCone.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

// ClassImp(TGeoVGSphere);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGSphere::TGeoVGSphere()
  : Base_t("", 0., 0., 0., 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoSph);
   fNz = 0;
   fNseg = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGSphere::TGeoVGSphere(Double_t rmin, Double_t rmax, Double_t theta1, Double_t theta2, Double_t phi1, Double_t phi2)
  : Base_t("", rmin, rmax, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad(), theta1 * TMath::DegToRad(), (theta2 - theta1) * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoSph);
   // SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGSphere::TGeoVGSphere(const char *name, Double_t rmin, Double_t rmax, Double_t theta1, Double_t theta2, Double_t phi1,
                       Double_t phi2)
  : Base_t(name, rmin, rmax, phi1 * TMath::DegToRad(), (phi2 - phi1) * TMath::DegToRad(), theta1 * TMath::DegToRad(), (theta2 - theta1) * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoSph);
   // SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
/// param[0] = Rmin
/// param[1] = Rmax
/// param[2] = theta1
/// param[3] = theta2
/// param[4] = phi1
/// param[5] = phi2

TGeoVGSphere::TGeoVGSphere(Double_t *param, Int_t /*nparam*/)
  : Base_t("", param[0], param[1], param[4]  * TMath::DegToRad(), (param[5] - param[4]) * TMath::DegToRad(), param[2]  * TMath::DegToRad(), (param[3] - param[2]) * TMath::DegToRad())
  {
   SetShapeBit(TGeoShape::kGeoSph);
   // SetDimensions(param, nparam);
   ComputeBBox();
   SetNumberOfDivisions(20);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGSphere::~TGeoVGSphere() {}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the sphere

void TGeoVGSphere::ComputeBBox()
{
   if (TGeoShape::IsSameWithinTolerance(TMath::Abs(GetDeltaThetaAngle() * TMath::RadToDeg()), 180)) {
      if (TGeoShape::IsSameWithinTolerance(TMath::Abs(GetDeltaPhiAngle() * TMath::RadToDeg()), 360)) {
         TGeoBBox::SetBoxDimensions(GetOuterRadius(), GetOuterRadius(), GetOuterRadius());
         memset(fOrigin, 0, 3 * sizeof(Double_t));
         return;
      }
   }
   Double_t st1 = TMath::Sin(GetStartThetaAngle());
   Double_t st2 = TMath::Sin(GetStartThetaAngle() + GetDeltaThetaAngle());
   Double_t r1min, r1max, r2min, r2max, rmin, rmax;
   r1min = TMath::Min(GetOuterRadius() * st1, GetOuterRadius() * st2);
   r1max = TMath::Max(GetOuterRadius() * st1, GetOuterRadius() * st2);
   r2min = TMath::Min(GetInnerRadius() * st1, GetInnerRadius() * st2);
   r2max = TMath::Max(GetInnerRadius() * st1, GetInnerRadius() * st2);
   if (((GetStartThetaAngle() * TMath::RadToDeg() <= 90) && ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() >= 90)) || (((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() <= 90) && (GetStartThetaAngle() * TMath::RadToDeg() >= 90))) {
      r1max = GetOuterRadius();
      r2max = GetInnerRadius();
   }
   rmin = TMath::Min(r1min, r2min);
   rmax = TMath::Max(r1max, r2max);

   Double_t xc[4];
   Double_t yc[4];
   xc[0] = rmax * TMath::Cos(GetStartPhiAngle());
   yc[0] = rmax * TMath::Sin(GetStartPhiAngle());
   xc[1] = rmax * TMath::Cos((GetStartPhiAngle() + GetDeltaPhiAngle()) );
   yc[1] = rmax * TMath::Sin((GetStartPhiAngle() + GetDeltaPhiAngle()) );
   xc[2] = rmin * TMath::Cos(GetStartPhiAngle());
   yc[2] = rmin * TMath::Sin(GetStartPhiAngle());
   xc[3] = rmin * TMath::Cos((GetStartPhiAngle() + GetDeltaPhiAngle()) );
   yc[3] = rmin * TMath::Sin((GetStartPhiAngle() + GetDeltaPhiAngle()) );

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];
   Double_t dp = (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg() - GetStartPhiAngle() * TMath::RadToDeg();
   if (dp < 0)
      dp += 360;
   Double_t ddp = -GetStartPhiAngle() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      xmax = rmax;
   ddp = 90 - GetStartPhiAngle() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      ymax = rmax;
   ddp = 180 - GetStartPhiAngle() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      xmin = -rmax;
   ddp = 270 - GetStartPhiAngle() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= dp)
      ymin = -rmax;
   xc[0] = GetOuterRadius() * TMath::Cos(GetStartThetaAngle());
   xc[1] = GetOuterRadius() * TMath::Cos((GetStartThetaAngle() + GetDeltaThetaAngle()));
   xc[2] = GetInnerRadius() * TMath::Cos(GetStartThetaAngle());
   xc[3] = GetInnerRadius() * TMath::Cos((GetStartThetaAngle() + GetDeltaThetaAngle()));
   Double_t zmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t zmax = xc[TMath::LocMax(4, &xc[0])];

   fOrigin[0] = (xmax + xmin) / 2;
   fOrigin[1] = (ymax + ymin) / 2;
   fOrigin[2] = (zmax + zmin) / 2;
   ;
   fDX = (xmax - xmin) / 2;
   fDY = (ymax - ymin) / 2;
   fDZ = (zmax - zmin) / 2;
}


////////////////////////////////////////////////////////////////////////////////

TGeoVolume *
TGeoVGSphere::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Int_t id;
   Double_t end = start + ndiv * step;
   switch (iaxis) {
   case 1: //---                R division
      finder = new TGeoPatternSphR(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         shape = new TGeoVGSphere(start + id * step, start + (id + 1) * step, GetStartThetaAngle() * TMath::RadToDeg(), (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg(), GetStartPhiAngle() * TMath::RadToDeg(), (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg());
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "R";
         voldiv->AddNodeOffset(vol, id, 0, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 2: //---                Phi division
      finder = new TGeoPatternSphPhi(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGSphere(GetInnerRadius(), GetOuterRadius(), GetStartThetaAngle() * TMath::RadToDeg(), (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg(), -step / 2, step / 2);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      opt = "Phi";
      for (id = 0; id < ndiv; id++) {
         voldiv->AddNodeOffset(vol, id, start + id * step + step / 2, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 3: //---                Theta division
      finder = new TGeoPatternSphTheta(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         shape = new TGeoVGSphere(GetInnerRadius(), GetOuterRadius(), start + id * step, start + (id + 1) * step, GetStartPhiAngle() * TMath::RadToDeg(), (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg());
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Theta";
         voldiv->AddNodeOffset(vol, id, 0, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   default: Error("Divide", "In shape %s wrong axis type for division", GetName()); return nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of axis IAXIS.

const char *TGeoVGSphere::GetAxisName(Int_t iaxis) const
{
   switch (iaxis) {
   case 1: return "R";
   case 2: return "PHI";
   case 3: return "THETA";
   default: return "UNDEFINED";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGSphere::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 1:
      xlo = GetInnerRadius();
      xhi = GetOuterRadius();
      dx = xhi - xlo;
      return dx;
   case 2:
      xlo = GetStartPhiAngle() * TMath::RadToDeg();
      xhi = (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg();
      dx = xhi - xlo;
      return dx;
   case 3:
      xlo = GetStartThetaAngle() * TMath::RadToDeg();
      xhi = (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg();
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGSphere::GetBoundingCylinder(Double_t *param) const
{
   Double_t smin = TMath::Sin(GetStartThetaAngle());
   Double_t smax = TMath::Sin((GetStartThetaAngle() + GetDeltaThetaAngle()));
   if (smin > smax) {
      Double_t a = smin;
      smin = smax;
      smax = a;
   }
   param[0] = GetInnerRadius() * smin; // Rmin
   param[0] *= param[0];
   if (((90. - GetStartThetaAngle() * TMath::RadToDeg()) * ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() - 90.)) >= 0)
      smax = 1.;
   param[1] = GetOuterRadius() * smax; // Rmax
   param[1] *= param[1];
   param[2] = (GetStartPhiAngle() * TMath::RadToDeg() < 0) ? (GetStartPhiAngle() * TMath::RadToDeg() + 360.) : GetStartPhiAngle() * TMath::RadToDeg(); // Phi1
   param[3] = (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg();
   if (TGeoShape::IsSameWithinTolerance(param[3] - param[2], 360)) { // Phi2
      param[2] = 0.;
      param[3] = 360.;
   }
   while (param[3] < param[2])
      param[3] += 360.;
}

////////////////////////////////////////////////////////////////////////////////
/// return start theta angle
Double_t TGeoVGSphere::GetTheta1() const
{
   return GetStartThetaAngle() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// return end theta angle
Double_t TGeoVGSphere::GetTheta2() const
{
   return (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// return start phi angle
Double_t TGeoVGSphere::GetPhi1() const
{
   return GetStartPhiAngle() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// return end phi angle
Double_t  TGeoVGSphere::GetPhi2() const
{
   return (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGSphere::InspectShape() const
{
   printf("*** Shape %s: TGeoVGSphere ***\n", GetName());
   printf("    Rmin = %11.5f\n", GetInnerRadius());
   printf("    Rmax = %11.5f\n", GetOuterRadius());
   printf("    Th1  = %11.5f\n", GetStartThetaAngle() * TMath::RadToDeg());
   printf("    Th2  = %11.5f\n", (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg());
   printf("    Ph1  = %11.5f\n", GetStartPhiAngle() * TMath::RadToDeg());
   printf("    Ph2  = %11.5f\n", (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGSphere::MakeBuffer3D() const
{
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
      full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg))
      ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;

   Int_t nbPnts = nlat * nlong + nup + ndown + ncenter;
   if (TestShapeBit(kGeoRSeg))
      nbPnts *= 2;

   Int_t nbSegs = nlat * fNseg + (nlat - 1 + nup + ndown) * nlong; // outer sphere
   if (TestShapeBit(kGeoRSeg))
      nbSegs *= 2; // inner sphere
   if (TestShapeBit(kGeoPhiSeg))
      nbSegs += 2 * nlat + nup + ndown; // 2 phi planes
   nbSegs += nlong * (2 - nup - ndown); // connecting cones

   Int_t nbPols = fNz * fNseg; // outer
   if (TestShapeBit(kGeoRSeg))
      nbPols *= 2; // inner
   if (TestShapeBit(kGeoPhiSeg))
      nbPols += 2 * fNz;                // 2 phi planes
   nbPols += (2 - nup - ndown) * fNseg; // connecting

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

void TGeoVGSphere::SetSegsAndPols(TBuffer3D &buff) const
{
   // Bool_t full = kTRUE;
   // if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg)) full = kFALSE;
   // Int_t ncenter = 1;
   // if (full || TestShapeBit(kGeoRSeg)) ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;

   // Int_t nbPnts = nlat*nlong+nup+ndown+ncenter;
   // if (TestShapeBit(kGeoRSeg)) nbPnts *= 2;

   // Int_t nbSegs = nlat*fNseg + (nlat-1+nup+ndown)*nlong; // outer sphere
   // if (TestShapeBit(kGeoRSeg)) nbSegs *= 2; // inner sphere
   // if (TestShapeBit(kGeoPhiSeg)) nbSegs += 2*nlat+nup+ndown; // 2 phi planes
   // nbSegs += nlong * (2-nup - ndown);  // connecting cones

   // Int_t nbPols = fNz*fNseg; // outer
   // if (TestShapeBit(kGeoRSeg)) nbPols *=2;  // inner
   // if (TestShapeBit(kGeoPhiSeg)) nbPols += 2*fNz; // 2 phi planes
   // nbPols += (2-nup-ndown)*fNseg; // connecting

   Int_t c = GetBasicColor();
   Int_t i, j;
   Int_t indx;
   indx = 0;
   // outside sphere
   // loop all segments on latitudes (except 0 and 180 degrees)
   // [0, nlat*fNseg)
   Int_t indpar = 0;
   for (i = 0; i < nlat; i++) {
      for (j = 0; j < fNseg; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = i * nlong + j;
         buff.fSegs[indx++] = i * nlong + (j + 1) % nlong;
      }
   }
   // loop all segments on longitudes
   // nlat*fNseg + [0, (nlat-1)*nlong)
   Int_t indlong = indpar + nlat * fNseg;
   for (i = 0; i < nlat - 1; i++) {
      for (j = 0; j < nlong; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = i * nlong + j;
         buff.fSegs[indx++] = (i + 1) * nlong + j;
      }
   }
   Int_t indup = indlong + (nlat - 1) * nlong;
   // extra longitudes on top
   // nlat*fNseg+(nlat-1)*nlong + [0, nlong)
   if (nup) {
      Int_t indpup = nlat * nlong;
      for (j = 0; j < nlong; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = j;
         buff.fSegs[indx++] = indpup;
      }
   }
   Int_t inddown = indup + nup * nlong;
   // extra longitudes on bottom
   // nlat*fNseg+(nlat+nup-1)*nlong + [0, nlong)
   if (ndown) {
      Int_t indpdown = nlat * nlong + nup;
      for (j = 0; j < nlong; j++) {
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = (nlat - 1) * nlong + j;
         buff.fSegs[indx++] = indpdown;
      }
   }
   Int_t indparin = inddown + ndown * nlong;
   Int_t indlongin = indparin;
   Int_t indupin = indparin;
   Int_t inddownin = indparin;
   Int_t indphi = indparin;
   // inner sphere
   Int_t indptin = nlat * nlong + nup + ndown;
   Int_t iptcenter = indptin;
   // nlat*fNseg+(nlat+nup+ndown-1)*nlong
   if (TestShapeBit(kGeoRSeg)) {
      indlongin = indparin + nlat * fNseg;
      indupin = indlongin + (nlat - 1) * nlong;
      inddownin = indupin + nup * nlong;
      // loop all segments on latitudes (except 0 and 180 degrees)
      // indsegin + [0, nlat*fNseg)
      for (i = 0; i < nlat; i++) {
         for (j = 0; j < fNseg; j++) {
            buff.fSegs[indx++] = c + 1;
            buff.fSegs[indx++] = indptin + i * nlong + j;
            buff.fSegs[indx++] = indptin + i * nlong + (j + 1) % nlong;
         }
      }
      // loop all segments on longitudes
      // indsegin + nlat*fNseg + [0, (nlat-1)*nlong)
      for (i = 0; i < nlat - 1; i++) {
         for (j = 0; j < nlong; j++) {
            buff.fSegs[indx++] = c + 1;
            buff.fSegs[indx++] = indptin + i * nlong + j;
            buff.fSegs[indx++] = indptin + (i + 1) * nlong + j;
         }
      }
      // extra longitudes on top
      // indsegin + nlat*fNseg+(nlat-1)*nlong + [0, nlong)
      if (nup) {
         Int_t indupltop = indptin + nlat * nlong;
         for (j = 0; j < nlong; j++) {
            buff.fSegs[indx++] = c + 1;
            buff.fSegs[indx++] = indptin + j;
            buff.fSegs[indx++] = indupltop;
         }
      }
      // extra longitudes on bottom
      // indsegin + nlat*fNseg+(nlat+nup-1)*nlong + [0, nlong)
      if (ndown) {
         Int_t indpdown = indptin + nlat * nlong + nup;
         for (j = 0; j < nlong; j++) {
            buff.fSegs[indx++] = c + 1;
            buff.fSegs[indx++] = indptin + (nlat - 1) * nlong + j;
            buff.fSegs[indx++] = indpdown;
         }
      }
      indphi = inddownin + ndown * nlong;
   }
   Int_t indtheta = indphi;
   // Segments on phi planes
   if (TestShapeBit(kGeoPhiSeg)) {
      indtheta += 2 * nlat + nup + ndown;
      for (j = 0; j < nlat; j++) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = j * nlong;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + j * nlong;
         else
            buff.fSegs[indx++] = iptcenter;
      }
      for (j = 0; j < nlat; j++) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = (j + 1) * nlong - 1;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + (j + 1) * nlong - 1;
         else
            buff.fSegs[indx++] = iptcenter;
      }
      if (nup) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = nlat * nlong;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + nlat * nlong;
         else
            buff.fSegs[indx++] = iptcenter;
      }
      if (ndown) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = nlat * nlong + nup;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + nlat * nlong + nup;
         else
            buff.fSegs[indx++] = iptcenter;
      }
   }
   // Segments on cones
   if (!nup) {
      for (j = 0; j < nlong; j++) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = j;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + j;
         else
            buff.fSegs[indx++] = iptcenter;
      }
   }
   if (!ndown) {
      for (j = 0; j < nlong; j++) {
         buff.fSegs[indx++] = c + 2;
         buff.fSegs[indx++] = (nlat - 1) * nlong + j;
         if (TestShapeBit(kGeoRSeg))
            buff.fSegs[indx++] = indptin + (nlat - 1) * nlong + j;
         else
            buff.fSegs[indx++] = iptcenter;
      }
   }

   indx = 0;
   // Fill polygons for outside sphere (except 0/180)
   for (i = 0; i < nlat - 1; i++) {
      for (j = 0; j < fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = indpar + i * fNseg + j;
         buff.fPols[indx++] = indlong + i * nlong + (j + 1) % nlong;
         buff.fPols[indx++] = indpar + (i + 1) * fNseg + j;
         buff.fPols[indx++] = indlong + i * nlong + j;
      }
   }
   // upper
   if (nup) {
      for (j = 0; j < fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 3;
         buff.fPols[indx++] = indup + j;
         buff.fPols[indx++] = indup + (j + 1) % nlong;
         buff.fPols[indx++] = indpar + j;
      }
   }
   // lower
   if (ndown) {
      for (j = 0; j < fNseg; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 3;
         buff.fPols[indx++] = inddown + j;
         buff.fPols[indx++] = indpar + (nlat - 1) * fNseg + j;
         buff.fPols[indx++] = inddown + (j + 1) % nlong;
      }
   }
   // Fill polygons for inside sphere (except 0/180)

   if (TestShapeBit(kGeoRSeg)) {
      for (i = 0; i < nlat - 1; i++) {
         for (j = 0; j < fNseg; j++) {
            buff.fPols[indx++] = c + 1;
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indparin + i * fNseg + j;
            buff.fPols[indx++] = indlongin + i * nlong + j;
            buff.fPols[indx++] = indparin + (i + 1) * fNseg + j;
            buff.fPols[indx++] = indlongin + i * nlong + (j + 1) % nlong;
         }
      }
      // upper
      if (nup) {
         for (j = 0; j < fNseg; j++) {
            buff.fPols[indx++] = c + 1;
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indupin + j;
            buff.fPols[indx++] = indparin + j;
            buff.fPols[indx++] = indupin + (j + 1) % nlong;
         }
      }
      // lower
      if (ndown) {
         for (j = 0; j < fNseg; j++) {
            buff.fPols[indx++] = c + 1;
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddownin + j;
            buff.fPols[indx++] = inddownin + (j + 1) % nlong;
            buff.fPols[indx++] = indparin + (nlat - 1) * fNseg + j;
         }
      }
   }
   // Polygons on phi planes
   if (TestShapeBit(kGeoPhiSeg)) {
      for (i = 0; i < nlat - 1; i++) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indlong + i * nlong;
            buff.fPols[indx++] = indphi + i + 1;
            buff.fPols[indx++] = indlongin + i * nlong;
            buff.fPols[indx++] = indphi + i;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indlong + i * nlong;
            buff.fPols[indx++] = indphi + i + 1;
            buff.fPols[indx++] = indphi + i;
         }
      }
      for (i = 0; i < nlat - 1; i++) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indlong + (i + 1) * nlong - 1;
            buff.fPols[indx++] = indphi + nlat + i;
            buff.fPols[indx++] = indlongin + (i + 1) * nlong - 1;
            buff.fPols[indx++] = indphi + nlat + i + 1;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indlong + (i + 1) * nlong - 1;
            buff.fPols[indx++] = indphi + nlat + i;
            buff.fPols[indx++] = indphi + nlat + i + 1;
         }
      }
      if (nup) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indup;
            buff.fPols[indx++] = indphi;
            buff.fPols[indx++] = indupin;
            buff.fPols[indx++] = indphi + 2 * nlat;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indup;
            buff.fPols[indx++] = indphi;
            buff.fPols[indx++] = indphi + 2 * nlat;
         }
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indup + nlong - 1;
            buff.fPols[indx++] = indphi + 2 * nlat;
            buff.fPols[indx++] = indupin + nlong - 1;
            buff.fPols[indx++] = indphi + nlat;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indup + nlong - 1;
            buff.fPols[indx++] = indphi + 2 * nlat;
            buff.fPols[indx++] = indphi + nlat;
         }
      }
      if (ndown) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = inddown;
            buff.fPols[indx++] = indphi + 2 * nlat + nup;
            buff.fPols[indx++] = inddownin;
            buff.fPols[indx++] = indphi + nlat - 1;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddown;
            buff.fPols[indx++] = indphi + 2 * nlat + nup;
            buff.fPols[indx++] = indphi + nlat - 1;
         }
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = inddown + nlong - 1;
            buff.fPols[indx++] = indphi + 2 * nlat - 1;
            buff.fPols[indx++] = inddownin + nlong - 1;
            buff.fPols[indx++] = indphi + 2 * nlat + nup;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = inddown + nlong - 1;
            buff.fPols[indx++] = indphi + 2 * nlat - 1;
            buff.fPols[indx++] = indphi + 2 * nlat + nup;
         }
      }
   }
   // Polygons on cones
   if (!nup) {
      for (j = 0; j < fNseg; j++) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indpar + j;
            buff.fPols[indx++] = indtheta + j;
            buff.fPols[indx++] = indparin + j;
            buff.fPols[indx++] = indtheta + (j + 1) % nlong;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indpar + j;
            buff.fPols[indx++] = indtheta + j;
            buff.fPols[indx++] = indtheta + (j + 1) % nlong;
         }
      }
   }
   if (!ndown) {
      for (j = 0; j < fNseg; j++) {
         buff.fPols[indx++] = c + 2;
         if (TestShapeBit(kGeoRSeg)) {
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indpar + (nlat - 1) * fNseg + j;
            buff.fPols[indx++] = indtheta + (1 - nup) * nlong + (j + 1) % nlong;
            buff.fPols[indx++] = indparin + (nlat - 1) * fNseg + j;
            buff.fPols[indx++] = indtheta + (1 - nup) * nlong + j;
         } else {
            buff.fPols[indx++] = 3;
            buff.fPols[indx++] = indpar + (nlat - 1) * fNseg + j;
            buff.fPols[indx++] = indtheta + (1 - nup) * nlong + (j + 1) % nlong;
            buff.fPols[indx++] = indtheta + (1 - nup) * nlong + j;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGSphere::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin   = " << GetInnerRadius() << ";" << std::endl;
   out << "   rmax   = " << GetOuterRadius() << ";" << std::endl;
   out << "   theta1 = " << GetStartThetaAngle() * TMath::RadToDeg() << ";" << std::endl;
   out << "   theta2 = " << (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() << ";" << std::endl;
   out << "   phi1   = " << GetStartPhiAngle() * TMath::RadToDeg() << ";" << std::endl;
   out << "   phi2   = " << (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGSphere(\"" << GetName()
       << "\",rmin,rmax,theta1, theta2,phi1,phi2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set spherical segment dimensions.

void TGeoVGSphere::SetSphDimensions(Double_t rmin, Double_t rmax, Double_t theta1, Double_t theta2, Double_t phi1,
                                  Double_t phi2)
{
   if (rmin >= rmax) {
      Error("SetDimensions", "invalid parameters rmin/rmax");
      return;
   }
   SetInnerRadius(rmin);
   SetOuterRadius(rmax);
   if (rmin > 0)
      SetShapeBit(kGeoRSeg);
   if (theta1 >= theta2 || theta1 < 0 || theta1 > 180 || theta2 > 180) {
      Error("SetDimensions", "invalid parameters theta1/theta2");
      return;
   }
   SetStartThetaAngle(theta1 * TMath::DegToRad());
   SetDeltaThetaAngle((theta2 - theta1) * TMath::DegToRad());
   if ((theta2 - theta1) < 180.)
      SetShapeBit(kGeoThetaSeg);
   auto newPhi1 = phi1;
   if (newPhi1 < 0)
      newPhi1 += 360.;
   SetStartPhiAngle(newPhi1 * TMath::DegToRad());
   auto newPhi2 = phi2;
   while (newPhi2 <= newPhi1)
      newPhi2 += 360.;
   SetDeltaPhiAngle((newPhi2 - newPhi1) * TMath::DegToRad());
   if (!TGeoShape::IsSameWithinTolerance(TMath::Abs(phi2 - phi1), 360))
   SetShapeBit(kGeoPhiSeg);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the spherical segment starting from a list of parameters.

void TGeoVGSphere::SetDimensions(Double_t *param, Int_t nparam)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t theta1 = 0;
   Double_t theta2 = 180.;
   Double_t phi1 = 0;
   Double_t phi2 = 360.;
   if (nparam > 2)
      theta1 = param[2];
   if (nparam > 3)
      theta2 = param[3];
   if (nparam > 4)
      phi1 = param[4];
   if (nparam > 5)
      phi2 = param[5];
   SetSphDimensions(rmin, rmax, theta1, theta2, phi1, phi2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the spherical segment starting from a list of parameters.
/// Only takes rmin and rmax

void TGeoVGSphere::SetDimensions(Double_t *param)
{
   SetDimensions(param, 2);
}

////////////////////////////////////////////////////////////////////////////////
/// Set the number of divisions of mesh circles keeping aspect ratio.

void TGeoVGSphere::SetNumberOfDivisions(Int_t p)
{
   fNseg = p;
   Double_t dphi = (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg() - GetStartPhiAngle() * TMath::RadToDeg();
   if (dphi < 0)
      dphi += 360;
   Double_t dtheta = TMath::Abs((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() - GetStartThetaAngle() * TMath::RadToDeg());
   fNz = Int_t(fNseg * dtheta / dphi) + 1;
   if (fNz < 2)
      fNz = 2;
}

////////////////////////////////////////////////////////////////////////////////
/// create sphere mesh points

void TGeoVGSphere::SetPoints(Double_t *points) const
{
   if (!points) {
      Error("SetPoints", "Input array is NULL");
      return;
   }
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
      full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg))
      ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t i, j;
   Double_t phi1 = GetStartPhiAngle();
   Double_t phi2 = (GetStartPhiAngle() + GetDeltaPhiAngle()) ;
   Double_t dphi = (phi2 - phi1) / fNseg;
   Double_t theta1 = GetStartThetaAngle();
   Double_t theta2 = (GetStartThetaAngle() + GetDeltaThetaAngle());
   Double_t dtheta = (theta2 - theta1) / fNz;
   Double_t z, zi, theta, phi, cphi, sphi;
   Int_t indx = 0;
   // FILL ALL POINTS ON OUTER SPHERE
   // (nlat * nlong) points
   // loop all latitudes except 0/180 degrees (nlat times)
   // ilat = [0,nlat]   jlong = [0,nlong]
   // Index(ilat, jlong) = 3*(ilat*nlat + jlong)
   for (i = 0; i < nlat; i++) {
      theta = theta1 + (nup + i) * dtheta;
      z = GetOuterRadius() * TMath::Cos(theta);
      zi = GetOuterRadius() * TMath::Sin(theta);
      // loop all different longitudes (nlong times)
      for (j = 0; j < nlong; j++) {
         phi = phi1 + j * dphi;
         cphi = TMath::Cos(phi);
         sphi = TMath::Sin(phi);
         points[indx++] = zi * cphi;
         points[indx++] = zi * sphi;
         points[indx++] = z;
      }
   }
   // upper/lower points (if they exist) for outer sphere
   if (nup) {
      // ind_up = 3*nlat*nlong
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = GetOuterRadius();
   }
   if (ndown) {
      // ind_down = 3*(nlat*nlong+nup)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -GetOuterRadius();
   }
   // do the same for inner sphere if it exist
   // Start_index = 3*(nlat*nlong + nup + ndown)
   if (TestShapeBit(kGeoRSeg)) {
      // Index(ilat, jlong) = start_index + 3*(ilat*nlat + jlong)
      for (i = 0; i < nlat; i++) {
         theta = theta1 + (nup + i) * dtheta;
         z = GetInnerRadius() * TMath::Cos(theta);
         zi = GetInnerRadius() * TMath::Sin(theta);
         // loop all different longitudes (nlong times)
         for (j = 0; j < nlong; j++) {
            phi = phi1 + j * dphi;
            cphi = TMath::Cos(phi);
            sphi = TMath::Sin(phi);
            points[indx++] = zi * cphi;
            points[indx++] = zi * sphi;
            points[indx++] = z;
         }
      }
      // upper/lower points (if they exist) for inner sphere
      if (nup) {
         // ind_up = start_index + 3*nlat*nlong
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = GetInnerRadius();
      }
      if (ndown) {
         // ind_down = start_index + 3*(nlat*nlong+nup)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -GetInnerRadius();
      }
   }
   // Add center of sphere if needed
   if (ncenter) {
      // ind_center = 6*(nlat*nlong + nup + ndown)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = 0.;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create sphere mesh points

void TGeoVGSphere::SetPoints(Float_t *points) const
{
   if (!points) {
      Error("SetPoints", "Input array is NULL");
      return;
   }
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
      full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg))
      ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t i, j;
   Double_t phi1 = GetStartPhiAngle();
   Double_t phi2 = (GetStartPhiAngle() + GetDeltaPhiAngle()) ;
   Double_t dphi = (phi2 - phi1) / fNseg;
   Double_t theta1 = GetStartThetaAngle();
   Double_t theta2 = (GetStartThetaAngle() + GetDeltaThetaAngle());
   Double_t dtheta = (theta2 - theta1) / fNz;
   Double_t z, zi, theta, phi, cphi, sphi;
   Int_t indx = 0;
   // FILL ALL POINTS ON OUTER SPHERE
   // (nlat * nlong) points
   // loop all latitudes except 0/180 degrees (nlat times)
   // ilat = [0,nlat]   jlong = [0,nlong]
   // Index(ilat, jlong) = 3*(ilat*nlat + jlong)
   for (i = 0; i < nlat; i++) {
      theta = theta1 + (nup + i) * dtheta;
      z = GetOuterRadius() * TMath::Cos(theta);
      zi = GetOuterRadius() * TMath::Sin(theta);
      // loop all different longitudes (nlong times)
      for (j = 0; j < nlong; j++) {
         phi = phi1 + j * dphi;
         cphi = TMath::Cos(phi);
         sphi = TMath::Sin(phi);
         points[indx++] = zi * cphi;
         points[indx++] = zi * sphi;
         points[indx++] = z;
      }
   }
   // upper/lower points (if they exist) for outer sphere
   if (nup) {
      // ind_up = 3*nlat*nlong
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = GetOuterRadius();
   }
   if (ndown) {
      // ind_down = 3*(nlat*nlong+nup)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -GetOuterRadius();
   }
   // do the same for inner sphere if it exist
   // Start_index = 3*(nlat*nlong + nup + ndown)
   if (TestShapeBit(kGeoRSeg)) {
      // Index(ilat, jlong) = start_index + 3*(ilat*nlat + jlong)
      for (i = 0; i < nlat; i++) {
         theta = theta1 + (nup + i) * dtheta;
         z = GetInnerRadius() * TMath::Cos(theta);
         zi = GetInnerRadius() * TMath::Sin(theta);
         // loop all different longitudes (nlong times)
         for (j = 0; j < nlong; j++) {
            phi = phi1 + j * dphi;
            cphi = TMath::Cos(phi);
            sphi = TMath::Sin(phi);
            points[indx++] = zi * cphi;
            points[indx++] = zi * sphi;
            points[indx++] = z;
         }
      }
      // upper/lower points (if they exist) for inner sphere
      if (nup) {
         // ind_up = start_index + 3*nlat*nlong
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = GetInnerRadius();
      }
      if (ndown) {
         // ind_down = start_index + 3*(nlat*nlong+nup)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -GetInnerRadius();
      }
   }
   // Add center of sphere if needed
   if (ncenter) {
      // ind_center = 6*(nlat*nlong + nup + ndown)
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = 0.;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGSphere::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   TGeoVGSphere *localThis = const_cast<TGeoVGSphere *>(this);
   localThis->SetNumberOfDivisions(gGeoManager->GetNsegments());
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
      full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg))
      ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;

   nvert = nlat * nlong + nup + ndown + ncenter;
   if (TestShapeBit(kGeoRSeg))
      nvert *= 2;

   nsegs = nlat * fNseg + (nlat - 1 + nup + ndown) * nlong; // outer sphere
   if (TestShapeBit(kGeoRSeg))
      nsegs *= 2; // inner sphere
   if (TestShapeBit(kGeoPhiSeg))
      nsegs += 2 * nlat + nup + ndown; // 2 phi planes
   nsegs += nlong * (2 - nup - ndown); // connecting cones

   npols = fNz * fNseg; // outer
   if (TestShapeBit(kGeoRSeg))
      npols *= 2; // inner
   if (TestShapeBit(kGeoPhiSeg))
      npols += 2 * fNz;                // 2 phi planes
   npols += (2 - nup - ndown) * fNseg; // connecting
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGSphere::GetNmeshVertices() const
{
   Bool_t full = kTRUE;
   if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
      full = kFALSE;
   Int_t ncenter = 1;
   if (full || TestShapeBit(kGeoRSeg))
      ncenter = 0;
   Int_t nup = (GetStartThetaAngle() * TMath::RadToDeg() > 0) ? 0 : 1;
   Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
   // number of different latitudes, excluding 0 and 180 degrees
   Int_t nlat = fNz + 1 - (nup + ndown);
   // number of different longitudes
   Int_t nlong = fNseg;
   if (TestShapeBit(kGeoPhiSeg))
      nlong++;
   // total number of points on mesh is:
   //    nlat*nlong + nup + ndown + ncenter;    // in case rmin=0
   //   2*(nlat*nlong + nup + ndown);           // in case rmin>0
   Int_t numPoints = 0;
   if (TestShapeBit(kGeoRSeg))
      numPoints = 2 * (nlat * nlong + nup + ndown);
   else
      numPoints = nlat * nlong + nup + ndown + ncenter;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
////// obsolete - to be removed

void TGeoVGSphere::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGSphere::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DSphere buffer;

   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusInner = GetInnerRadius();
      buffer.fRadiusOuter = GetOuterRadius();
      buffer.fThetaMin = GetStartThetaAngle() * TMath::RadToDeg();
      buffer.fThetaMax = (GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg();
      buffer.fPhiMin = GetStartPhiAngle() * TMath::RadToDeg();
      buffer.fPhiMax = (GetStartPhiAngle() + GetDeltaPhiAngle()) * TMath::RadToDeg();
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      // We want FillBuffer to be const
      TGeoVGSphere *localThis = const_cast<TGeoVGSphere *>(this);
      localThis->SetNumberOfDivisions(gGeoManager->GetNsegments());

      Bool_t full = kTRUE;
      if (TestShapeBit(kGeoThetaSeg) || TestShapeBit(kGeoPhiSeg))
         full = kFALSE;
      Int_t ncenter = 1;
      if (full || TestShapeBit(kGeoRSeg))
         ncenter = 0;
      Int_t nup = (GetStartThetaAngle()* TMath::RadToDeg() > 0) ? 0 : 1;
      Int_t ndown = ((GetStartThetaAngle() + GetDeltaThetaAngle()) * TMath::RadToDeg() < 180) ? 0 : 1;
      // number of different latitudes, excluding 0 and 180 degrees
      Int_t nlat = fNz + 1 - (nup + ndown);
      // number of different longitudes
      Int_t nlong = fNseg;
      if (TestShapeBit(kGeoPhiSeg))
         nlong++;

      Int_t nbPnts = nlat * nlong + nup + ndown + ncenter;
      if (TestShapeBit(kGeoRSeg))
         nbPnts *= 2;

      Int_t nbSegs = nlat * fNseg + (nlat - 1 + nup + ndown) * nlong; // outer sphere
      if (TestShapeBit(kGeoRSeg))
         nbSegs *= 2; // inner sphere
      if (TestShapeBit(kGeoPhiSeg))
         nbSegs += 2 * nlat + nup + ndown; // 2 phi planes
      nbSegs += nlong * (2 - nup - ndown); // connecting cones

      Int_t nbPols = fNz * fNseg; // outer
      if (TestShapeBit(kGeoRSeg))
         nbPols *= 2; // inner
      if (TestShapeBit(kGeoPhiSeg))
         nbPols += 2 * fNz;                // 2 phi planes
      nbPols += (2 - nup - ndown) * fNseg; // connecting

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

#endif
