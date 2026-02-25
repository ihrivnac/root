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

// ClassImp(TGeoVGGtra);

////////////////////////////////////////////////////////////////////////////////
/// Default ctor

TGeoVGGtra::TGeoVGGtra()
  : Base_t("")
{
   fTwistAngle = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoVGGtra::TGeoVGGtra(Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1, Double_t bl1, Double_t tl1,
                   Double_t alpha1, Double_t h2, Double_t bl2, Double_t tl2, Double_t alpha2)
   : Base_t("")
{
   // From TGeoTrap base
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

   // From TGeoGtra
   fTwistAngle = twist;
   Double_t x, y;
   Double_t th = theta * TMath::DegToRad();
   Double_t ph = phi * TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -dz * TMath::Sin(th) * TMath::Cos(ph);
   Double_t yc = -dz * TMath::Sin(th) * TMath::Sin(ph);

   Int_t i;
   for (i = 0; i < 4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i = 4; i < 8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(0.5 * twist * TMath::DegToRad()) + yc;
   }
   // VG shape
   Double_t verticesx[8], verticesy[8];
   for (Int_t i = 0; i < 8; i++) {
      verticesx[i] = fXY[i][0];
      verticesy[i] = fXY[i][1];
   }
   Initialize(verticesx, verticesy, dz);

   ComputeTwist();
   if ((dz < 0) || (h1 < 0) || (bl1 < 0) || (tl1 < 0) || (h2 < 0) || (bl2 < 0) || (tl2 < 0))
      SetShapeBit(kGeoRunTimeShape);
   else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor providing the name of the shape.

TGeoVGGtra::TGeoVGGtra(const char *name, Double_t dz, Double_t theta, Double_t phi, Double_t twist, Double_t h1,
                   Double_t bl1, Double_t tl1, Double_t alpha1, Double_t h2, Double_t bl2, Double_t tl2,
                   Double_t alpha2)
   : Base_t(name)
{
   // From TGeoTrap base
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

   // From TGeoGtra
   fTwistAngle = twist;
   Double_t x, y;
   Double_t th = theta * TMath::DegToRad();
   Double_t ph = phi * TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -dz * TMath::Sin(th) * TMath::Cos(ph);
   Double_t yc = -dz * TMath::Sin(th) * TMath::Sin(ph);

   Int_t i;

   for (i = 0; i < 4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i = 4; i < 8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(0.5 * twist * TMath::DegToRad()) + yc;
   }
   // VG shape
   Double_t verticesx[8], verticesy[8];
   for (Int_t i = 0; i < 8; i++) {
      verticesx[i] = fXY[i][0];
      verticesy[i] = fXY[i][1];
   }
   Initialize(verticesx, verticesy, dz);
   ComputeTwist();
   if ((dz < 0) || (h1 < 0) || (bl1 < 0) || (tl1 < 0) || (h2 < 0) || (bl2 < 0) || (tl2 < 0))
      SetShapeBit(kGeoRunTimeShape);
   else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVGGtra::~TGeoVGGtra()
{
    if (fTwist)
    delete[] fTwist;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy twist values from source array
/// From TGeoArb8

void TGeoVGGtra::CopyTwist(Double_t *twist)
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
/// Computes bounding box for an Arb8 shape.
/// From TGeoArb8

void TGeoVGGtra::ComputeBBox()
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

   Double_t dx, dy, dz;
   Double_t origin[3];
   dx = 0.5 * (xmax - xmin);
   dy = 0.5 * (ymax - ymin);
   dz = fDz;
   origin[0] = 0.5 * (xmax + xmin);
   origin[1] = 0.5 * (ymax + ymin);
   origin[2] = 0;
   fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);
   SetShapeBit(kGeoClosedShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes tangents of twist angles (angles between projections on XY plane
/// of corresponding -dz +dz edges). Computes also if the vertices are defined
/// clockwise or anti-clockwise.
/// From TGeoArb8

void TGeoVGGtra::ComputeTwist()
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
/// Divide this shape along one axis.

TGeoVolume *TGeoVGGtra::Divide(TGeoVolume *voldiv, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Division of a general trapezoid not implemented");
   return voldiv;
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGGtra::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz     = " << fDz << ";" << std::endl;
   out << "   theta  = " << fTheta << ";" << std::endl;
   out << "   phi    = " << fPhi << ";" << std::endl;
   out << "   twist  = " << fTwistAngle << ";" << std::endl;
   out << "   h1     = " << fH1 << ";" << std::endl;
   out << "   bl1    = " << fBl1 << ";" << std::endl;
   out << "   tl1    = " << fTl1 << ";" << std::endl;
   out << "   alpha1 = " << fAlpha1 << ";" << std::endl;
   out << "   h2     = " << fH2 << ";" << std::endl;
   out << "   bl2    = " << fBl2 << ";" << std::endl;
   out << "   tl2    = " << fTl2 << ";" << std::endl;
   out << "   alpha2 = " << fAlpha2 << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGGtra(\"" << GetName()
       << "\", dz,theta,phi,twist,h1,bl1,tl1,alpha1,h2,bl2,tl2,alpha2);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGGtra::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (mother->IsRunTimeShape()) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t dz, h1, bl1, tl1, h2, bl2, tl2;
   if (fDz < 0)
      dz = ((TGeoVGGtra *)mother)->GetDz();
   else
      dz = fDz;
   if (fH1 < 0)
      h1 = ((TGeoVGGtra *)mother)->GetH1();
   else
      h1 = fH1;
   if (fH2 < 0)
      h2 = ((TGeoVGGtra *)mother)->GetH2();
   else
      h2 = fH2;
   if (fBl1 < 0)
      bl1 = ((TGeoVGGtra *)mother)->GetBl1();
   else
      bl1 = fBl1;
   if (fBl2 < 0)
      bl2 = ((TGeoVGGtra *)mother)->GetBl2();
   else
      bl2 = fBl2;
   if (fTl1 < 0)
      tl1 = ((TGeoVGGtra *)mother)->GetTl1();
   else
      tl1 = fTl1;
   if (fTl2 < 0)
      tl2 = ((TGeoVGGtra *)mother)->GetTl2();
   else
      tl2 = fTl2;
   return (new TGeoVGGtra(dz, fTheta, fPhi, fTwistAngle, h1, bl1, tl1, fAlpha1, h2, bl2, tl2, fAlpha2));
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
///  - param[11] = twist

void TGeoVGGtra::SetDimensions(Double_t *param)
{
   // From TGeoTrap::SetDimensions(param);
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

   // From TGeoGtra::SetDimensions(param);
   fTwistAngle = param[11];
   Double_t x, y;
   Double_t twist = fTwistAngle;
   Double_t th = fTheta * TMath::DegToRad();
   Double_t ph = fPhi * TMath::DegToRad();
   // Coordinates of the center of the bottom face
   Double_t xc = -fDz * TMath::Sin(th) * TMath::Cos(ph);
   Double_t yc = -fDz * TMath::Sin(th) * TMath::Sin(ph);

   Int_t i;
   for (i = 0; i < 4; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(-0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(-0.5 * twist * TMath::DegToRad()) + yc;
   }
   // Coordinates of the center of the top face
   xc = -xc;
   yc = -yc;
   for (i = 4; i < 8; i++) {
      x = fXY[i][0] - xc;
      y = fXY[i][1] - yc;
      fXY[i][0] =
         x * TMath::Cos(0.5 * twist * TMath::DegToRad()) + y * TMath::Sin(0.5 * twist * TMath::DegToRad()) + xc;
      fXY[i][1] =
         -x * TMath::Sin(0.5 * twist * TMath::DegToRad()) + y * TMath::Cos(0.5 * twist * TMath::DegToRad()) + yc;
   }
   // VG shape
   Double_t verticesx[8], verticesy[8];
   for (Int_t i = 0; i < 8; i++) {
      verticesx[i] = fXY[i][0];
      verticesy[i] = fXY[i][1];
   }
   Initialize(verticesx, verticesy, fDz);
   ComputeTwist();
   if ((fDz < 0) || (fH1 < 0) || (fBl1 < 0) || (fTl1 < 0) || (fH2 < 0) || (fBl2 < 0) || (fTl2 < 0))
      SetShapeBit(kGeoRunTimeShape);
   else
      ComputeBBox();
}

#endif
