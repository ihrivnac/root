// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoVGArb8.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>
#include "TBuffer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TMath.h"

// ClassImp(TGeoVGTrap);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor

TGeoVGTrap::TGeoVGTrap()
  : Base_t("", 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
{
   fDz = 0;
   fTheta = 0;
   fPhi = 0;
   fH1 = fH2 = fBl1 = fBl2 = fTl1 = fTl2 = fAlpha1 = fAlpha2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing just a range in Z, theta and phi.
/// To DO: check if it makes sense

TGeoVGTrap::TGeoVGTrap(Double_t dz, Double_t theta, Double_t phi)
  : Base_t("", dz, theta * TMath::DegToRad(), phi  * TMath::DegToRad(), 0., 0., 0., 0., 0., 0., 0., 0.)
{
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = fH2 = fBl1 = fBl2 = fTl1 = fTl2 = fAlpha1 = fAlpha2 = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Normal constructor.

TGeoVGTrap::TGeoVGTrap(Double_t dz, Double_t theta, Double_t phi, Double_t h1, Double_t bl1, Double_t tl1, Double_t alpha1,
                   Double_t h2, Double_t bl2, Double_t tl2, Double_t alpha2)
   : Base_t("", dz, theta  * TMath::DegToRad(), phi * TMath::DegToRad(), 
            h1, bl1, tl1, alpha1 * TMath::DegToRad(), h2, bl2, tl2, alpha2 * TMath::DegToRad())
{
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = h1;
   fH2 = h2;
   fBl1 = bl1;
   fBl2 = bl2;
   fTl1 = tl1;
   fTl2 = tl2;
   fAlpha1 = alpha1;
   fAlpha2 = alpha2;
   Double_t tx = TMath::Tan(theta * TMath::DegToRad()) * TMath::Cos(phi * TMath::DegToRad());
   Double_t ty = TMath::Tan(theta * TMath::DegToRad()) * TMath::Sin(phi * TMath::DegToRad());
   Double_t ta1 = TMath::Tan(alpha1 * TMath::DegToRad());
   Double_t ta2 = TMath::Tan(alpha2 * TMath::DegToRad());
   fXY[0][0] = -dz * tx - h1 * ta1 - bl1;
   fXY[0][1] = -dz * ty - h1;
   fXY[1][0] = -dz * tx + h1 * ta1 - tl1;
   fXY[1][1] = -dz * ty + h1;
   fXY[2][0] = -dz * tx + h1 * ta1 + tl1;
   fXY[2][1] = -dz * ty + h1;
   fXY[3][0] = -dz * tx - h1 * ta1 + bl1;
   fXY[3][1] = -dz * ty - h1;
   fXY[4][0] = dz * tx - h2 * ta2 - bl2;
   fXY[4][1] = dz * ty - h2;
   fXY[5][0] = dz * tx + h2 * ta2 - tl2;
   fXY[5][1] = dz * ty + h2;
   fXY[6][0] = dz * tx + h2 * ta2 + tl2;
   fXY[6][1] = dz * ty + h2;
   fXY[7][0] = dz * tx - h2 * ta2 + bl2;
   fXY[7][1] = dz * ty - h2;
   ComputeTwist();
   if ((dz < 0) || (h1 < 0) || (bl1 < 0) || (tl1 < 0) || (h2 < 0) || (bl2 < 0) || (tl2 < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with name.

TGeoVGTrap::TGeoVGTrap(const char *name, Double_t dz, Double_t theta, Double_t phi, Double_t h1, Double_t bl1, Double_t tl1,
                   Double_t alpha1, Double_t h2, Double_t bl2, Double_t tl2, Double_t alpha2)
   : Base_t(name, dz, theta  * TMath::DegToRad(), phi * TMath::DegToRad(), 
            h1, bl1, tl1, alpha1 * TMath::DegToRad(), h2, bl2, tl2, alpha2 * TMath::DegToRad())
{
   fDz = dz;
   fTheta = theta;
   fPhi = phi;
   fH1 = h1;
   fH2 = h2;
   fBl1 = bl1;
   fBl2 = bl2;
   fTl1 = tl1;
   fTl2 = tl2;
   fAlpha1 = alpha1;
   fAlpha2 = alpha2;
   for (Int_t i = 0; i < 8; i++) {
      fXY[i][0] = 0.0;
      fXY[i][1] = 0.0;
   }
   Double_t tx = TMath::Tan(theta * TMath::DegToRad()) * TMath::Cos(phi * TMath::DegToRad());
   Double_t ty = TMath::Tan(theta * TMath::DegToRad()) * TMath::Sin(phi * TMath::DegToRad());
   Double_t ta1 = TMath::Tan(alpha1 * TMath::DegToRad());
   Double_t ta2 = TMath::Tan(alpha2 * TMath::DegToRad());
   fXY[0][0] = -dz * tx - h1 * ta1 - bl1;
   fXY[0][1] = -dz * ty - h1;
   fXY[1][0] = -dz * tx + h1 * ta1 - tl1;
   fXY[1][1] = -dz * ty + h1;
   fXY[2][0] = -dz * tx + h1 * ta1 + tl1;
   fXY[2][1] = -dz * ty + h1;
   fXY[3][0] = -dz * tx - h1 * ta1 + bl1;
   fXY[3][1] = -dz * ty - h1;
   fXY[4][0] = dz * tx - h2 * ta2 - bl2;
   fXY[4][1] = dz * ty - h2;
   fXY[5][0] = dz * tx + h2 * ta2 - tl2;
   fXY[5][1] = dz * ty + h2;
   fXY[6][0] = dz * tx + h2 * ta2 + tl2;
   fXY[6][1] = dz * ty + h2;
   fXY[7][0] = dz * tx - h2 * ta2 + bl2;
   fXY[7][1] = dz * ty - h2;
   ComputeTwist();
   if ((dz < 0) || (h1 < 0) || (bl1 < 0) || (tl1 < 0) || (h2 < 0) || (bl2 < 0) || (tl2 < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVGTrap::~TGeoVGTrap() {}

////////////////////////////////////////////////////////////////////////////////
/// Computes bounding box.
/// From TGeoArb8

void TGeoVGTrap::ComputeBBox()
{
   Double_t xmin, xmax, ymin, ymax;
   xmin = xmax = fXY[0][0];
   ymin = ymax = fXY[0][1];

   for (Int_t i = 1; i < 8; i++) {
      if (xmin > fXY[i][0])
         xmin = fXY[i][0];
      if (xmax < fXY[i][0])
         xmax = fXY[i][0];
      if (ymin > fXY[i][1])
         ymin = fXY[i][1];
      if (ymax < fXY[i][1])
         ymax = fXY[i][1];
   }
   fDX = 0.5 * (xmax - xmin);
   fDY = 0.5 * (ymax - ymin);
   fDZ = fDz;
   fOrigin[0] = 0.5 * (xmax + xmin);
   fOrigin[1] = 0.5 * (ymax + ymin);
   fOrigin[2] = 0;
   SetShapeBit(kGeoClosedShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes tangents of twist angles (angles between projections on XY plane
/// of corresponding -dz +dz edges). Computes also if the vertices are defined
/// clockwise or anti-clockwise.
/// From TGeoArb8

void TGeoVGTrap::ComputeTwist()
{
   Double_t twist[4];
   Bool_t twisted = kFALSE;
   Double_t dx1, dy1, dx2, dy2;
   Bool_t singleBottom = kTRUE;
   Bool_t singleTop = kTRUE;
   Int_t i;
   for (i = 0; i < 4; i++) {
      dx1 = fXY[(i + 1) % 4][0] - fXY[i][0];
      dy1 = fXY[(i + 1) % 4][1] - fXY[i][1];
      if (TMath::Abs(dx1) < TGeoShape::Tolerance() && TMath::Abs(dy1) < TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      singleBottom = kFALSE;
      dx2 = fXY[4 + (i + 1) % 4][0] - fXY[4 + i][0];
      dy2 = fXY[4 + (i + 1) % 4][1] - fXY[4 + i][1];
      if (TMath::Abs(dx2) < TGeoShape::Tolerance() && TMath::Abs(dy2) < TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      singleTop = kFALSE;
      twist[i] = dy1 * dx2 - dx1 * dy2;
      if (TMath::Abs(twist[i]) < TGeoShape::Tolerance()) {
         twist[i] = 0;
         continue;
      }
      twist[i] = TMath::Sign(1., twist[i]);
      twisted = kTRUE;
   }

   CopyTwist(twisted ? twist : nullptr);

   if (singleBottom) {
      for (i = 0; i < 4; i++) {
         fXY[i][0] += 1.E-8 * fXY[i + 4][0];
         fXY[i][1] += 1.E-8 * fXY[i + 4][1];
      }
   }
   if (singleTop) {
      for (i = 0; i < 4; i++) {
         fXY[i + 4][0] += 1.E-8 * fXY[i][0];
         fXY[i + 4][1] += 1.E-8 * fXY[i][1];
      }
   }
   Double_t sum1 = 0.;
   Double_t sum2 = 0.;
   Int_t j;
   for (i = 0; i < 4; i++) {
      j = (i + 1) % 4;
      sum1 += fXY[i][0] * fXY[j][1] - fXY[j][0] * fXY[i][1];
      sum2 += fXY[i + 4][0] * fXY[j + 4][1] - fXY[j + 4][0] * fXY[i + 4][1];
   }
   if (sum1 * sum2 < -TGeoShape::Tolerance()) {
      Fatal("ComputeTwist", "Shape %s type Arb8: Lower/upper faces defined with opposite clockwise", GetName());
      return;
   }
   if (sum1 > TGeoShape::Tolerance()) {
      Error("ComputeTwist", "Shape %s type Arb8: Vertices must be defined clockwise in XY planes. Re-ordering...",
            GetName());
      Double_t xtemp, ytemp;
      xtemp = fXY[1][0];
      ytemp = fXY[1][1];
      fXY[1][0] = fXY[3][0];
      fXY[1][1] = fXY[3][1];
      fXY[3][0] = xtemp;
      fXY[3][1] = ytemp;
      xtemp = fXY[5][0];
      ytemp = fXY[5][1];
      fXY[5][0] = fXY[7][0];
      fXY[5][1] = fXY[7][1];
      fXY[7][0] = xtemp;
      fXY[7][1] = ytemp;
   }
   // Check for illegal crossings.
   Bool_t illegal_cross = kFALSE;
   illegal_cross =
      TGeoShape::IsSegCrossing(fXY[0][0], fXY[0][1], fXY[1][0], fXY[1][1], fXY[2][0], fXY[2][1], fXY[3][0], fXY[3][1]);
   if (!illegal_cross)
      illegal_cross = TGeoShape::IsSegCrossing(fXY[4][0], fXY[4][1], fXY[5][0], fXY[5][1], fXY[6][0], fXY[6][1],
                                               fXY[7][0], fXY[7][1]);
   if (illegal_cross) {
      Error("ComputeTwist", "Shape %s type Arb8: Malformed polygon with crossing opposite segments", GetName());
      InspectShape();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Copy twist values from source array
/// From TGeoArb8

void TGeoVGTrap::CopyTwist(Double_t *twist)
{
   if (twist) {
      if (!fTwist)
         fTwist = new Double_t[4];
      memcpy(fTwist, twist, 4 * sizeof(Double_t));
   } else if (fTwist) {
      delete[] fTwist;
      fTwist = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Computes intersection points between plane at zpl and non-horizontal edges.
/// From TGeoArb8

void TGeoVGTrap::SetPlaneVertices(Double_t zpl, Double_t *vertices) const
/// From TGeoArb8
{
   Double_t cf = 0.5 * (fDz - zpl) / fDz;
   for (Int_t i = 0; i < 4; i++) {
      vertices[2 * i] = fXY[i + 4][0] + cf * (fXY[i][0] - fXY[i + 4][0]);
      vertices[2 * i + 1] = fXY[i + 4][1] + cf * (fXY[i][1] - fXY[i + 4][1]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this trapezoid shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Only Z divisions
/// are supported. For Z divisions just return the pointer to the volume to be
/// divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *
TGeoVGTrap::Divide(TGeoVolume */*voldiv*/, const char */*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/, Double_t /*start*/, Double_t /*step*/)
{
/*
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   if (iaxis != 3) {
      Error("Divide", "cannot divide trapezoids on other axis than Z");
      return nullptr;
   }
   Double_t end = start + ndiv * step;
   Double_t points_lo[8];
   Double_t points_hi[8];
   finder = new TGeoPatternTrapZ(voldiv, ndiv, start, end);
   voldiv->SetFinder(finder);
   finder->SetDivIndex(voldiv->GetNdaughters());
   opt = "Z";
   vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
   Double_t txz = ((TGeoPatternTrapZ *)finder)->GetTxz();
   Double_t tyz = ((TGeoPatternTrapZ *)finder)->GetTyz();
   Double_t zmin, zmax, ox, oy, oz;
   for (Int_t idiv = 0; idiv < ndiv; idiv++) {
      zmin = start + idiv * step;
      zmax = start + (idiv + 1) * step;
      oz = start + idiv * step + step / 2;
      ox = oz * txz;
      oy = oz * tyz;
      SetPlaneVertices(zmin, &points_lo[0]);
      SetPlaneVertices(zmax, &points_hi[0]);
      shape = new TGeoVGTrap(step / 2, fTheta, fPhi);
      for (Int_t vert1 = 0; vert1 < 4; vert1++)
         ((TGeoArb8 *)shape)->SetVertex(vert1, points_lo[2 * vert1] - ox, points_lo[2 * vert1 + 1] - oy);
      for (Int_t vert2 = 0; vert2 < 4; vert2++)
         ((TGeoArb8 *)shape)->SetVertex(vert2 + 4, points_hi[2 * vert2] - ox, points_hi[2 * vert2 + 1] - oy);
      vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
      vmulti->AddVolume(vol);
      voldiv->AddNodeOffset(vol, idiv, oz, opt.Data());
      ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
   }
   return vmulti;
*/
   Error("Divide", "Cannot divide VG trapezoid");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these have to be computed
/// in order to fit the mother.

TGeoShape *TGeoVGTrap::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (mother->IsRunTimeShape()) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t dz, h1, bl1, tl1, h2, bl2, tl2;
   if (fDz < 0)
      dz = ((TGeoVGTrap *)mother)->GetDz();
   else
      dz = fDz;

   if (fH1 < 0)
      h1 = ((TGeoVGTrap *)mother)->GetH1();
   else
      h1 = fH1;

   if (fH2 < 0)
      h2 = ((TGeoVGTrap *)mother)->GetH2();
   else
      h2 = fH2;

   if (fBl1 < 0)
      bl1 = ((TGeoVGTrap *)mother)->GetBl1();
   else
      bl1 = fBl1;

   if (fBl2 < 0)
      bl2 = ((TGeoVGTrap *)mother)->GetBl2();
   else
      bl2 = fBl2;

   if (fTl1 < 0)
      tl1 = ((TGeoVGTrap *)mother)->GetTl1();
   else
      tl1 = fTl1;

   if (fTl2 < 0)
      tl2 = ((TGeoVGTrap *)mother)->GetTl2();
   else
      tl2 = fTl2;

   return (new TGeoVGTrap(dz, fTheta, fPhi, h1, bl1, tl1, fAlpha1, h2, bl2, tl2, fAlpha2));
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGTrap::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz     = " << fDz << ";" << std::endl;
   out << "   theta  = " << fTheta << ";" << std::endl;
   out << "   phi    = " << fPhi << ";" << std::endl;
   out << "   h1     = " << fH1 << ";" << std::endl;
   out << "   bl1    = " << fBl1 << ";" << std::endl;
   out << "   tl1    = " << fTl1 << ";" << std::endl;
   out << "   alpha1 = " << fAlpha1 << ";" << std::endl;
   out << "   h2     = " << fH2 << ";" << std::endl;
   out << "   bl2    = " << fBl2 << ";" << std::endl;
   out << "   tl2    = " << fTl2 << ";" << std::endl;
   out << "   alpha2 = " << fAlpha2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGTrap(\"" << GetName()
       << "\", dz,theta,phi,h1,bl1,tl1,alpha1,h2,bl2,tl2,alpha2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set all arb8 params in one step.
///  - param[0] = dz
///  - param[1] = theta
///  - param[2] = phi
///  - param[3] = h1
///  - param[4] = bl1
///  - param[5] = tl1
///  - param[6] = alpha1
///  - param[7] = h2
///  - param[8] = bl2
///  - param[9] = tl2
///  - param[10] = alpha2

void TGeoVGTrap::SetDimensions(Double_t *param)
{
   fDz = param[0];
   fTheta = param[1];
   fPhi = param[2];
   fH1 = param[3];
   fH2 = param[7];
   fBl1 = param[4];
   fBl2 = param[8];
   fTl1 = param[5];
   fTl2 = param[9];
   fAlpha1 = param[6];
   fAlpha2 = param[10];
   Double_t tx = TMath::Tan(fTheta * TMath::DegToRad()) * TMath::Cos(fPhi * TMath::DegToRad());
   Double_t ty = TMath::Tan(fTheta * TMath::DegToRad()) * TMath::Sin(fPhi * TMath::DegToRad());
   Double_t ta1 = TMath::Tan(fAlpha1 * TMath::DegToRad());
   Double_t ta2 = TMath::Tan(fAlpha2 * TMath::DegToRad());
   fXY[0][0] = -fDz * tx - fH1 * ta1 - fBl1;
   fXY[0][1] = -fDz * ty - fH1;
   fXY[1][0] = -fDz * tx + fH1 * ta1 - fTl1;
   fXY[1][1] = -fDz * ty + fH1;
   fXY[2][0] = -fDz * tx + fH1 * ta1 + fTl1;
   fXY[2][1] = -fDz * ty + fH1;
   fXY[3][0] = -fDz * tx - fH1 * ta1 + fBl1;
   fXY[3][1] = -fDz * ty - fH1;
   fXY[4][0] = fDz * tx - fH2 * ta2 - fBl2;
   fXY[4][1] = fDz * ty - fH2;
   fXY[5][0] = fDz * tx + fH2 * ta2 - fTl2;
   fXY[5][1] = fDz * ty + fH2;
   fXY[6][0] = fDz * tx + fH2 * ta2 + fTl2;
   fXY[6][1] = fDz * ty + fH2;
   fXY[7][0] = fDz * tx - fH2 * ta2 + fBl2;
   fXY[7][1] = fDz * ty - fH2;
   ComputeTwist();
   if ((fDz < 0) || (fH1 < 0) || (fBl1 < 0) || (fTl1 < 0) || (fH2 < 0) || (fBl2 < 0) || (fTl2 < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   } else
      ComputeBBox();
}

#endif
