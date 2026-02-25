// @(#)root/geom:$Id$
// Author: Andrei Gheata   28/07/03

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoTorus.h"
#include "TGeoVGTorus.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoVGTorus.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGTorus);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGTorus::TGeoVGTorus()
   : Base_t("", 0., 0., 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoTorus);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor without name.

TGeoVGTorus::TGeoVGTorus(Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
   : Base_t("", rmin, rmax, r, phi1 * TMath::DegToRad(), dphi * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetTorusDimensions(r, rmin, rmax, phi1, dphi);
   if ((rmin < 0) || (rmax < 0))
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with name.

TGeoVGTorus::TGeoVGTorus(const char *name, Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
   : Base_t(name, rmin, rmax, r, phi1 * TMath::DegToRad(), dphi * TMath::DegToRad())
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetTorusDimensions(r, rmin, rmax, phi1, dphi);
   if ((rmin < 0) || (rmax < 0))
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor based on an array of parameters.
///  - param[0] = R
///  - param[1] = Rmin
///  - param[2] = Rmax
///  - param[3] = Phi1
///  - param[4] = Dphi

TGeoVGTorus::TGeoVGTorus(Double_t *param)
   : Base_t("", 0., 0., 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoTorus);
   SetDimensions(param);
   if (param[1] < 0 || param[2] < 0)
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box of the torus.

void TGeoVGTorus::ComputeBBox()
{
   Double_t dx, dy, dz;
   Double_t origin[3];
   dz = rmax();   
   if (TGeoShape::IsSameWithinTolerance(dphi() * TMath::RadToDeg(), 360)) {
      dx = dy = rtor() + rmax();
      origin[0] = 0.;
      origin[1] = 0.;
      origin[2] = 0;
      fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);
      return;
   }
   Double_t xc[4];
   Double_t yc[4];
   xc[0] = (rtor() + rmax()) * TMath::Cos(sphi());
   yc[0] = (rtor() + rmax()) * TMath::Sin(sphi());
   xc[1] = (rtor() + rmax()) * TMath::Cos(sphi() + dphi());
   yc[1] = (rtor() + rmax()) * TMath::Sin(sphi() + dphi());
   xc[2] = (rtor() - rmax()) * TMath::Cos(sphi());
   yc[2] = (rtor() - rmax()) * TMath::Sin(sphi());
   xc[3] = (rtor() - rmax()) * TMath::Cos(sphi() + dphi());
   yc[3] = (rtor() - rmax()) * TMath::Sin(sphi() + dphi());

   Double_t xmin = xc[TMath::LocMin(4, &xc[0])];
   Double_t xmax = xc[TMath::LocMax(4, &xc[0])];
   Double_t ymin = yc[TMath::LocMin(4, &yc[0])];
   Double_t ymax = yc[TMath::LocMax(4, &yc[0])];
   Double_t ddp = -sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp <= (dphi() * TMath::RadToDeg()))
      xmax = rtor() + rmax();
   ddp = 90 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= (dphi() * TMath::RadToDeg()))
      ymax = rtor() + rmax();
   ddp = 180 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= (dphi() * TMath::RadToDeg()))
      xmin = -(rtor() + rmax());
   ddp = 270 - sphi() * TMath::RadToDeg();
   if (ddp < 0)
      ddp += 360;
   if (ddp > 360)
      ddp -= 360;
   if (ddp <= (dphi() * TMath::RadToDeg()))
      ymin = -(rtor() + rmax());

   dx = (xmax - xmin) / 2;
   dy = (ymax - ymin) / 2;
   origin[0] = (xmax + xmin) / 2;
   origin[1] = (ymax + ymin) / 2;
   origin[2] = 0;
   fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);
}


////////////////////////////////////////////////////////////////////////////////
/// Computes distance to axis of the torus from point pt + t*dir;

Double_t TGeoVGTorus::Daxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i = 0; i < 3; i++)
      p[i] = pt[i] + t * dir[i];
   Double_t rxy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
   return TMath::Sqrt((rxy - rtor()) * (rxy - rtor()) + p[2] * p[2]);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes derivative w.r.t. t of the distance to axis of the torus from point pt + t*dir;

Double_t TGeoVGTorus::DDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i = 0; i < 3; i++)
      p[i] = pt[i] + t * dir[i];
   Double_t rxy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
   if (rxy < 1E-4)
      return ((p[2] * dir[2] - rtor() * TMath::Sqrt(dir[0] * dir[0] + dir[1] * dir[1])) /
              TMath::Sqrt(rtor() * rtor() + p[2] * p[2]));
   Double_t d = TMath::Sqrt((rxy - rtor()) * (rxy - rtor()) + p[2] * p[2]);
   if (TGeoShape::IsSameWithinTolerance(d, 0))
      return 0.;
   Double_t dd = (p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2] - (p[0] * dir[0] + p[1] * dir[1]) * rtor() / rxy) / d;
   return dd;
}

////////////////////////////////////////////////////////////////////////////////
/// Second derivative of distance to torus axis w.r.t t.

Double_t TGeoVGTorus::DDDaxis(const Double_t *pt, const Double_t *dir, Double_t t) const
{
   Double_t p[3];
   for (Int_t i = 0; i < 3; i++)
      p[i] = pt[i] + t * dir[i];
   Double_t rxy = TMath::Sqrt(p[0] * p[0] + p[1] * p[1]);
   if (rxy < 1E-6)
      return 0;
   Double_t daxis = TMath::Sqrt((rxy - rtor()) * (rxy - rtor()) + p[2] * p[2]);
   if (TGeoShape::IsSameWithinTolerance(daxis, 0))
      return 0;
   Double_t ddaxis =
      (p[0] * dir[0] + p[1] * dir[1] + p[2] * dir[2] - (p[0] * dir[0] + p[1] * dir[1]) * rtor() / rxy) / daxis;
   Double_t dddaxis = 1 - ddaxis * ddaxis - (1 - dir[2] * dir[2]) * rtor() / rxy +
                      rtor() * (p[0] * dir[0] + p[1] * dir[1]) * (p[0] * dir[0] + p[1] * dir[1]) / (rxy * rxy * rxy);
   dddaxis /= daxis;
   return dddaxis;
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this torus shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step.

TGeoVolume *TGeoVGTorus::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                              Double_t /*start*/, Double_t /*step*/)
{
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns name of axis IAXIS.

const char *TGeoVGTorus::GetAxisName(Int_t iaxis) const
{
   switch (iaxis) {
   case 1: return "R";
   case 2: return "PHI";
   case 3: return "Z";
   default: return "UNDEFINED";
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGTorus::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
      dx = dphi() * TMath::RadToDeg();
      return dx;
   case 3: dx = 0; return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2, dZ

void TGeoVGTorus::GetBoundingCylinder(Double_t *param) const
{
   param[0] = (rtor() - rmax());  // Rmin
   param[1] = (rtor() + rmax());  // Rmax
   param[2] = sphi() * TMath::RadToDeg();         // Phi1
   param[3] = (sphi() + dphi()) * TMath::RadToDeg(); // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// Create a shape fitting the mother.

TGeoShape *TGeoVGTorus::GetMakeRuntimeShape(TGeoShape * /*mother*/, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   Error("GetMakeRuntimeShape", "parametrized toruses not supported");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGTorus::InspectShape() const
{
   printf("*** Shape %s: TGeoVGTorus ***\n", GetName());
   printf("    R    = %11.5f\n", rtor());
   printf("    Rmin = %11.5f\n", rmin());
   printf("    Rmax = %11.5f\n", rmax());
   printf("    Phi1 = %11.5f\n", sphi() * TMath::RadToDeg());
   printf("    Dphi = %11.5f\n", dphi() * TMath::RadToDeg());
   printf(" Bounding box:\n");
   fBoundingBox.InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGTorus::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t nbPnts = n * (n - 1);
   Bool_t hasrmin = (GetRmin() > 0) ? kTRUE : kFALSE;
   Bool_t hasphi = (GetDphi() < 360) ? kTRUE : kFALSE;
   if (hasrmin)
      nbPnts *= 2;
   else if (hasphi)
      nbPnts += 2;

   Int_t nbSegs = (2 * n - 1) * (n - 1);
   Int_t nbPols = (n - 1) * (n - 1);
   if (hasrmin) {
      nbSegs += (2 * n - 1) * (n - 1);
      nbPols += (n - 1) * (n - 1);
   }
   if (hasphi) {
      nbSegs += 2 * (n - 1);
      nbPols += 2 * (n - 1);
   }

   TBuffer3D *buff =
      new TBuffer3D(TBuffer3DTypes::kGeneric, nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols);
   if (buff) {
      SetPoints(buff->fPnts);
      SetSegsAndPols(*buff);
   }

   return buff;
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi1

Double_t TGeoVGTorus::GetPhi1() const
{
   return sphi() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Return phi2

Double_t TGeoVGTorus::GetDphi() const
{
   return dphi() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons.

void TGeoVGTorus::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t i, j;
   Int_t n = gGeoManager->GetNsegments() + 1;
   // Int_t nbPnts = n*(n-1);
   Int_t indx, indp, startcap = 0;
   Bool_t hasrmin = (GetRmin() > 0) ? kTRUE : kFALSE;
   Bool_t hasphi = (GetDphi() < 360) ? kTRUE : kFALSE;
   // if (hasrmin) nbPnts *= 2;
   // else if (hasphi) nbPnts += 2;
   Int_t c = GetBasicColor();

   indp = n * (n - 1); // start index for points on inner surface
   memset(buff.fSegs, 0, buff.NbSegs() * 3 * sizeof(Int_t));

   // outer surface phi circles = n*(n-1) -> [0, n*(n-1) -1]
   // connect point j with point j+1 on same row
   indx = 0;
   for (i = 0; i < n; i++) {        // rows [0,n-1]
      for (j = 0; j < n - 1; j++) { // points on a row [0, n-2]
         buff.fSegs[indx + (i * (n - 1) + j) * 3] = c;
         buff.fSegs[indx + (i * (n - 1) + j) * 3 + 1] = i * (n - 1) + j;                   // j on row i
         buff.fSegs[indx + (i * (n - 1) + j) * 3 + 2] = i * (n - 1) + ((j + 1) % (n - 1)); // j+1 on row i
      }
   }
   indx += 3 * n * (n - 1);
   // outer surface generators = (n-1)*(n-1) -> [n*(n-1), (2*n-1)*(n-1) -1]
   // connect point j on row i with point j on row i+1
   for (i = 0; i < n - 1; i++) {    // rows [0, n-2]
      for (j = 0; j < n - 1; j++) { // points on a row [0, n-2]
         buff.fSegs[indx + (i * (n - 1) + j) * 3] = c;
         buff.fSegs[indx + (i * (n - 1) + j) * 3 + 1] = i * (n - 1) + j;       // j on row i
         buff.fSegs[indx + (i * (n - 1) + j) * 3 + 2] = (i + 1) * (n - 1) + j; // j on row i+1
      }
   }
   indx += 3 * (n - 1) * (n - 1);
   startcap = (2 * n - 1) * (n - 1);

   if (hasrmin) {
      // inner surface phi circles = n*(n-1) -> [(2*n-1)*(n-1), (3*n-1)*(n-1) -1]
      // connect point j with point j+1 on same row
      for (i = 0; i < n; i++) {                                                    // rows [0, n-1]
         for (j = 0; j < n - 1; j++) {                                             // points on a row [0, n-2]
            buff.fSegs[indx + (i * (n - 1) + j) * 3] = c;                          // lighter color
            buff.fSegs[indx + (i * (n - 1) + j) * 3 + 1] = indp + i * (n - 1) + j; // j on row i
            buff.fSegs[indx + (i * (n - 1) + j) * 3 + 2] = indp + i * (n - 1) + ((j + 1) % (n - 1)); // j+1 on row i
         }
      }
      indx += 3 * n * (n - 1);
      // inner surface generators = (n-1)*n -> [(3*n-1)*(n-1), (4*n-2)*(n-1) -1]
      // connect point j on row i with point j on row i+1
      for (i = 0; i < n - 1; i++) {                                                      // rows [0, n-2]
         for (j = 0; j < n - 1; j++) {                                                   // points on a row [0, n-2]
            buff.fSegs[indx + (i * (n - 1) + j) * 3] = c;                                // lighter color
            buff.fSegs[indx + (i * (n - 1) + j) * 3 + 1] = indp + i * (n - 1) + j;       // j on row i
            buff.fSegs[indx + (i * (n - 1) + j) * 3 + 2] = indp + (i + 1) * (n - 1) + j; // j on row i+1
         }
      }
      indx += 3 * (n - 1) * (n - 1);
      startcap = (4 * n - 2) * (n - 1);
   }

   if (hasphi) {
      if (hasrmin) {
         // endcaps = 2*(n-1) -> [(4*n-2)*(n-1), 4*n*(n-1)-1]
         i = 0;
         for (j = 0; j < n - 1; j++) {
            buff.fSegs[indx + j * 3] = c + 1;
            buff.fSegs[indx + j * 3 + 1] = (n - 1) * i + j;        // outer j on row 0
            buff.fSegs[indx + j * 3 + 2] = indp + (n - 1) * i + j; // inner j on row 0
         }
         indx += 3 * (n - 1);
         i = n - 1;
         for (j = 0; j < n - 1; j++) {
            buff.fSegs[indx + j * 3] = c + 1;
            buff.fSegs[indx + j * 3 + 1] = (n - 1) * i + j;        // outer j on row n-1
            buff.fSegs[indx + j * 3 + 2] = indp + (n - 1) * i + j; // inner j on row n-1
         }
         indx += 3 * (n - 1);
      } else {
         i = 0;
         for (j = 0; j < n - 1; j++) {
            buff.fSegs[indx + j * 3] = c + 1;
            buff.fSegs[indx + j * 3 + 1] = (n - 1) * i + j; // outer j on row 0
            buff.fSegs[indx + j * 3 + 2] = n * (n - 1);     // center of first endcap
         }
         indx += 3 * (n - 1);
         i = n - 1;
         for (j = 0; j < n - 1; j++) {
            buff.fSegs[indx + j * 3] = c + 1;
            buff.fSegs[indx + j * 3 + 1] = (n - 1) * i + j; // outer j on row n-1
            buff.fSegs[indx + j * 3 + 2] = n * (n - 1) + 1; // center of second endcap
         }
         indx += 3 * (n - 1);
      }
   }

   indx = 0;
   memset(buff.fPols, 0, buff.NbPols() * 6 * sizeof(Int_t));

   // outer surface = (n-1)*(n-1) -> [0, (n-1)*(n-1)-1]
   // normal pointing out
   for (i = 0; i < n - 1; i++) {
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = n * (n - 1) + (n - 1) * i + ((j + 1) % (n - 1)); // generator j+1 on outer row i
         buff.fPols[indx++] = (n - 1) * (i + 1) + j;                           // seg j on outer row i+1
         buff.fPols[indx++] = n * (n - 1) + (n - 1) * i + j;                   // generator j on outer row i
         buff.fPols[indx++] = (n - 1) * i + j;                                 // seg j on outer row i
      }
   }
   if (hasrmin) {
      indp = (2 * n - 1) * (n - 1); // start index of inner segments
      // inner surface = (n-1)*(n-1) -> [(n-1)*(n-1), 2*(n-1)*(n-1)-1]
      // normal pointing out
      for (i = 0; i < n - 1; i++) {
         for (j = 0; j < n - 1; j++) {
            buff.fPols[indx++] = c;
            buff.fPols[indx++] = 4;
            buff.fPols[indx++] = indp + n * (n - 1) + (n - 1) * i + j;                   // generator j on inner row i
            buff.fPols[indx++] = indp + (n - 1) * (i + 1) + j;                           // seg j on inner row i+1
            buff.fPols[indx++] = indp + n * (n - 1) + (n - 1) * i + ((j + 1) % (n - 1)); // generator j+1 on inner r>
            buff.fPols[indx++] = indp + (n - 1) * i + j;                                 // seg j on inner row i
         }
      }
   }
   if (hasphi) {
      // endcaps = 2*(n-1) -> [2*(n-1)*(n-1), 2*n*(n-1)-1]
      i = 0; // row 0
      Int_t np = (hasrmin) ? 4 : 3;
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c + 1;
         buff.fPols[indx++] = np;
         buff.fPols[indx++] = j;            // seg j on outer row 0  a
         buff.fPols[indx++] = startcap + j; // endcap j on row 0  d
         if (hasrmin)
            buff.fPols[indx++] = indp + j;                    // seg j on inner row 0  c
         buff.fPols[indx++] = startcap + ((j + 1) % (n - 1)); // endcap j+1 on row 0  b
      }

      i = n - 1; // row n-1
      for (j = 0; j < n - 1; j++) {
         buff.fPols[indx++] = c + 1;
         buff.fPols[indx++] = np;
         buff.fPols[indx++] = (n - 1) * i + j;                          // seg j on outer row n-1 a
         buff.fPols[indx++] = startcap + (n - 1) + ((j + 1) % (n - 1)); // endcap j+1 on row n-1 d
         if (hasrmin)
            buff.fPols[indx++] = indp + (n - 1) * i + j; // seg j on inner row n-1 c
         buff.fPols[indx++] = startcap + (n - 1) + j;    // endcap j on row n-1 b
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGTorus::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   r    = " << rtor() << ";" << std::endl;
   out << "   rmin = " << rmin() << ";" << std::endl;
   out << "   rmax = " << rmax() << ";" << std::endl;
   out << "   phi1 = " << sphi() * TMath::RadToDeg() << ";" << std::endl;
   out << "   dphi = " << dphi() * TMath::RadToDeg() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGTorus(\"" << GetName() << "\",r,rmin,rmax,phi1,dphi);"
       << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set torus dimensions.

void TGeoVGTorus::SetTorusDimensions(Double_t r, Double_t rmin, Double_t rmax, Double_t phi1, Double_t dphi)
{
   SetRMin(rmin);
   SetRMax(rmax);
   SetRTor(r);
   auto phiStart = phi1;
   if (phiStart < 0)
      phiStart += 360.;
   SetSPhi(phiStart* TMath::DegToRad());
   SetDPhi(dphi * TMath::DegToRad());
}

////////////////////////////////////////////////////////////////////////////////
/// Set torus dimensions starting from a list.

void TGeoVGTorus::SetDimensions(Double_t *param)
{
   SetTorusDimensions(param[0], param[1], param[2], param[3], param[4]);
}

////////////////////////////////////////////////////////////////////////////////
/// Create torus mesh points

void TGeoVGTorus::SetPoints(Double_t *points) const
{
   if (!points)
      return;
   Int_t n = gGeoManager->GetNsegments() + 1;
   Double_t phin, phout;
   Double_t dpin = 360. / (n - 1);
   Double_t dpout = dphi() * TMath::RadToDeg() / (n - 1);
   Double_t co, so, ci, si;
   Bool_t havermin = (rmin() < TGeoShape::Tolerance()) ? kFALSE : kTRUE;
   Int_t i, j;
   Int_t indx = 0;
   // loop outer mesh -> n*(n-1) points [0, 3*n*(n-1)-1]
   for (i = 0; i < n; i++) {
      phout = (sphi() * TMath::RadToDeg() + i * dpout) * TMath::DegToRad();
      co = TMath::Cos(phout);
      so = TMath::Sin(phout);
      for (j = 0; j < n - 1; j++) {
         phin = j * dpin * TMath::DegToRad();
         ci = TMath::Cos(phin);
         si = TMath::Sin(phin);
         points[indx++] = (rtor() + rmax() * ci) * co;
         points[indx++] = (rtor() + rmax() * ci) * so;
         points[indx++] = rmax() * si;
      }
   }

   if (havermin) {
      // loop inner mesh -> n*(n-1) points [3*n*(n-1), 6*n*(n-1)]
      for (i = 0; i < n; i++) {
         phout = (sphi() * TMath::RadToDeg() + i * dpout) * TMath::DegToRad();
         co = TMath::Cos(phout);
         so = TMath::Sin(phout);
         for (j = 0; j < n - 1; j++) {
            phin = j * dpin * TMath::DegToRad();
            ci = TMath::Cos(phin);
            si = TMath::Sin(phin);
            points[indx++] = (rtor() + rmin() * ci) * co;
            points[indx++] = (rtor() + rmin() * ci) * so;
            points[indx++] = rmin() * si;
         }
      }
   } else {
      if ((dphi() * TMath::RadToDeg()) < 360.) {
         // just add extra 2 points on the centers of the 2 phi cuts [3*n*n, 3*n*n+1]
         points[indx++] = rtor() * TMath::Cos(sphi());
         points[indx++] = rtor() * TMath::Sin(sphi());
         points[indx++] = 0;
         points[indx++] = rtor() * TMath::Cos(sphi() + dphi());
         points[indx++] = rtor() * TMath::Sin(sphi() + dphi());
         points[indx++] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Create torus mesh points

void TGeoVGTorus::SetPoints(Float_t *points) const
{
   if (!points)
      return;
   Int_t n = gGeoManager->GetNsegments() + 1;
   Double_t phin, phout;
   Double_t dpin = 360. / (n - 1);
   Double_t dpout = (dphi() * TMath::RadToDeg()) / (n - 1);
   Double_t co, so, ci, si;
   Bool_t havermin = (rmin() < TGeoShape::Tolerance()) ? kFALSE : kTRUE;
   Int_t i, j;
   Int_t indx = 0;
   // loop outer mesh -> n*(n-1) points [0, 3*n*(n-1)-1]
   // plane i = 0, n-1  point j = 0, n-1  ipoint = n*i + j
   for (i = 0; i < n; i++) {
      phout = (sphi() * TMath::RadToDeg() + i * dpout) * TMath::DegToRad();
      co = TMath::Cos(phout);
      so = TMath::Sin(phout);
      for (j = 0; j < n - 1; j++) {
         phin = j * dpin * TMath::DegToRad();
         ci = TMath::Cos(phin);
         si = TMath::Sin(phin);
         points[indx++] = (rtor() + rmax() * ci) * co;
         points[indx++] = (rtor() + rmax() * ci) * so;
         points[indx++] = rmax() * si;
      }
   }

   if (havermin) {
      // loop inner mesh -> n*(n-1) points [3*n*(n-1), 6*n*(n-1)]
      // plane i = 0, n-1  point j = 0, n-1  ipoint = n*n + n*i + j
      for (i = 0; i < n; i++) {
         phout = (sphi() * TMath::RadToDeg() + i * dpout) * TMath::DegToRad();
         co = TMath::Cos(phout);
         so = TMath::Sin(phout);
         for (j = 0; j < n - 1; j++) {
            phin = j * dpin * TMath::DegToRad();
            ci = TMath::Cos(phin);
            si = TMath::Sin(phin);
            points[indx++] = (rtor() + rmin() * ci) * co;
            points[indx++] = (rtor() + rmin() * ci) * so;
            points[indx++] = rmin() * si;
         }
      }
   } else {
      if ((dphi() * TMath::RadToDeg()) < 360.) {
         // just add extra 2 points on the centers of the 2 phi cuts [n*n, n*n+1]
         // ip1 = n*(n-1) + 0;
         // ip2 = n*(n-1) + 1
         points[indx++] = rtor() * TMath::Cos(sphi());
         points[indx++] = rtor() * TMath::Sin(sphi());
         points[indx++] = 0;
         points[indx++] = rtor() * TMath::Cos(sphi() + dphi());
         points[indx++] = rtor() * TMath::Sin(sphi() + dphi());
         points[indx++] = 0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGTorus::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t numPoints = n * (n - 1);
   if (rmin() > TGeoShape::Tolerance())
      numPoints *= 2;
   else if ((dphi() * TMath::RadToDeg()) < 360.)
      numPoints += 2;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGTorus::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Find real solutions of the cubic equation : x^3 + a*x^2 + b*x + c = 0
/// Input: a,b,c
/// Output: x[3] real solutions
/// Returns number of real solutions (1 or 3)

Int_t TGeoVGTorus::SolveCubic(Double_t a, Double_t b, Double_t c, Double_t *x) const
{
   const Double_t ott = 1. / 3.;
   const Double_t sq3 = TMath::Sqrt(3.);
   Int_t ireal = 1;
   Double_t p = b - a * a * ott;
   Double_t q = c - a * b * ott + 2. * a * a * a * ott * ott * ott;
   Double_t delta = 4 * p * p * p + 27 * q * q;
   //   Double_t y1r, y1i, y2r, y2i;
   Double_t t, u;
   if (delta >= 0) {
      delta = TMath::Sqrt(delta);
      t = (-3 * q * sq3 + delta) / (6 * sq3);
      u = (3 * q * sq3 + delta) / (6 * sq3);
      x[0] = TMath::Sign(1., t) * TMath::Power(TMath::Abs(t), ott) -
             TMath::Sign(1., u) * TMath::Power(TMath::Abs(u), ott) - a * ott;
   } else {
      delta = TMath::Sqrt(-delta);
      t = -0.5 * q;
      u = delta / (6 * sq3);
      x[0] = 2. * TMath::Power(t * t + u * u, 0.5 * ott) * TMath::Cos(ott * TMath::ATan2(u, t));
      x[0] -= a * ott;
   }

   t = x[0] * x[0] + a * x[0] + b;
   u = a + x[0];
   delta = u * u - 4. * t;
   if (delta >= 0) {
      ireal = 3;
      delta = TMath::Sqrt(delta);
      x[1] = 0.5 * (-u - delta);
      x[2] = 0.5 * (-u + delta);
   }
   return ireal;
}

////////////////////////////////////////////////////////////////////////////////
/// Find real solutions of the quartic equation : x^4 + a*x^3 + b*x^2 + c*x + d = 0
/// Input: a,b,c,d
/// Output: x[4] - real solutions
/// Returns number of real solutions (0 to 3)

Int_t TGeoVGTorus::SolveQuartic(Double_t a, Double_t b, Double_t c, Double_t d, Double_t *x) const
{
   Double_t e = b - 3. * a * a / 8.;
   Double_t f = c + a * a * a / 8. - 0.5 * a * b;
   Double_t g = d - 3. * a * a * a * a / 256. + a * a * b / 16. - a * c / 4.;
   Double_t xx[4];
   Int_t ind[4];
   Double_t delta;
   Double_t h = 0;
   Int_t ireal = 0;
   Int_t i;
   if (TGeoShape::IsSameWithinTolerance(f, 0)) {
      delta = e * e - 4. * g;
      if (delta < 0)
         return 0;
      delta = TMath::Sqrt(delta);
      h = 0.5 * (-e - delta);
      if (h >= 0) {
         h = TMath::Sqrt(h);
         x[ireal++] = -h - 0.25 * a;
         x[ireal++] = h - 0.25 * a;
      }
      h = 0.5 * (-e + delta);
      if (h >= 0) {
         h = TMath::Sqrt(h);
         x[ireal++] = -h - 0.25 * a;
         x[ireal++] = h - 0.25 * a;
      }
      if (ireal > 0) {
         TMath::Sort(ireal, x, ind, kFALSE);
         for (i = 0; i < ireal; i++)
            xx[i] = x[ind[i]];
         memcpy(x, xx, ireal * sizeof(Double_t));
      }
      return ireal;
   }

   if (TGeoShape::IsSameWithinTolerance(g, 0)) {
      x[ireal++] = -0.25 * a;
      ind[0] = SolveCubic(0, e, f, xx);
      for (i = 0; i < ind[0]; i++)
         x[ireal++] = xx[i] - 0.25 * a;
      if (ireal > 0) {
         TMath::Sort(ireal, x, ind, kFALSE);
         for (i = 0; i < ireal; i++)
            xx[i] = x[ind[i]];
         memcpy(x, xx, ireal * sizeof(Double_t));
      }
      return ireal;
   }

   ireal = SolveCubic(2. * e, e * e - 4. * g, -f * f, xx);
   if (ireal == 1) {
      if (xx[0] <= 0)
         return 0;
      h = TMath::Sqrt(xx[0]);
   } else {
      // 3 real solutions of the cubic
      for (i = 0; i < 3; i++) {
         h = xx[i];
         if (h >= 0)
            break;
      }
      if (h <= 0)
         return 0;
      h = TMath::Sqrt(h);
   }
   Double_t j = 0.5 * (e + h * h - f / h);
   ireal = 0;
   delta = h * h - 4. * j;
   if (delta >= 0) {
      delta = TMath::Sqrt(delta);
      x[ireal++] = 0.5 * (-h - delta) - 0.25 * a;
      x[ireal++] = 0.5 * (-h + delta) - 0.25 * a;
   }
   delta = h * h - 4. * g / j;
   if (delta >= 0) {
      delta = TMath::Sqrt(delta);
      x[ireal++] = 0.5 * (h - delta) - 0.25 * a;
      x[ireal++] = 0.5 * (h + delta) - 0.25 * a;
   }
   if (ireal > 0) {
      TMath::Sort(ireal, x, ind, kFALSE);
      for (i = 0; i < ireal; i++)
         xx[i] = x[ind[i]];
      memcpy(x, xx, ireal * sizeof(Double_t));
   }
   return ireal;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns distance to the surface or the torus (fR,r) from a point, along
/// a direction. Point is close enough to the boundary so that the distance
/// to the torus is decreasing while moving along the given direction.

Double_t TGeoVGTorus::ToBoundary(const Double_t *pt, const Double_t *dir, Double_t r, Bool_t in) const
{
   // Compute coefficients of the quartic
   Double_t tol = TGeoShape::Tolerance();
   Double_t r0sq = pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2];
   Double_t rdotn = pt[0] * dir[0] + pt[1] * dir[1] + pt[2] * dir[2];
   Double_t rsumsq = rtor() * rtor() + r * r;
   Double_t a = 4. * rdotn;
   Double_t b = 2. * (r0sq + 2. * rdotn * rdotn - rsumsq + 2. * rtor() * rtor() * dir[2] * dir[2]);
   Double_t c = 4. * (r0sq * rdotn - rsumsq * rdotn + 2. * rtor() * rtor() * pt[2] * dir[2]);
   Double_t d = r0sq * r0sq - 2. * r0sq * rsumsq + 4. * rtor() * rtor() * pt[2] * pt[2] + (rtor() * rtor() - r * r) * (rtor() * rtor() - r * r);

   Double_t x[4], y[4];
   Int_t nsol = 0;

   if (TMath::Abs(dir[2]) < 1E-3 && TMath::Abs(pt[2]) < 0.1 * r) {
      Double_t r0 = rtor() - TMath::Sqrt((r - pt[2]) * (r + pt[2]));
      Double_t b0 = (pt[0] * dir[0] + pt[1] * dir[1]) / (dir[0] * dir[0] + dir[1] * dir[1]);
      Double_t c0 = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) / (dir[0] * dir[0] + dir[1] * dir[1]);
      Double_t delta = b0 * b0 - c0;
      if (delta > 0) {
         y[nsol] = -b0 - TMath::Sqrt(delta);
         if (y[nsol] > -tol)
            nsol++;
         y[nsol] = -b0 + TMath::Sqrt(delta);
         if (y[nsol] > -tol)
            nsol++;
      }
      r0 = rtor() + TMath::Sqrt((r - pt[2]) * (r + pt[2]));
      c0 = (pt[0] * pt[0] + (pt[1] - r0) * (pt[1] + r0)) / (dir[0] * dir[0] + dir[1] * dir[1]);
      delta = b0 * b0 - c0;
      if (delta > 0) {
         y[nsol] = -b0 - TMath::Sqrt(delta);
         if (y[nsol] > -tol)
            nsol++;
         y[nsol] = -b0 + TMath::Sqrt(delta);
         if (y[nsol] > -tol)
            nsol++;
      }
      if (nsol) {
         // Sort solutions
         Int_t ind[4];
         TMath::Sort(nsol, y, ind, kFALSE);
         for (Int_t j = 0; j < nsol; j++)
            x[j] = y[ind[j]];
      }
   } else {
      nsol = SolveQuartic(a, b, c, d, x);
   }
   if (!nsol)
      return TGeoShape::Big();
   // look for first positive solution
   Double_t phi, ndotd;
   Double_t r0[3], norm[3];
   Bool_t inner = (TMath::Abs(r - rmin()) < TGeoShape::Tolerance()) ? kTRUE : kFALSE;
   for (Int_t i = 0; i < nsol; i++) {
      if (x[i] < -10)
         continue;
      phi = TMath::ATan2(pt[1] + x[i] * dir[1], pt[0] + x[i] * dir[0]);
      r0[0] = rtor() * TMath::Cos(phi);
      r0[1] = rtor() * TMath::Sin(phi);
      r0[2] = 0;
      for (Int_t ipt = 0; ipt < 3; ipt++)
         norm[ipt] = pt[ipt] + x[i] * dir[ipt] - r0[ipt];
      ndotd = norm[0] * dir[0] + norm[1] * dir[1] + norm[2] * dir[2];
      if (inner ^ in) {
         if (ndotd < 0)
            continue;
      } else {
         if (ndotd > 0)
            continue;
      }
      Double_t s = x[i];
      Double_t eps = TGeoShape::Big();
      Double_t delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
      Double_t eps0 = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
      while (TMath::Abs(eps) > TGeoShape::Tolerance()) {
         if (TMath::Abs(eps0) > 100)
            break;
         s += eps0;
         if (TMath::Abs(s + eps0) < TGeoShape::Tolerance())
            break;
         delta = s * s * s * s + a * s * s * s + b * s * s + c * s + d;
         eps = -delta / (4. * s * s * s + 3. * a * s * s + 2. * b * s + c);
         if (TMath::Abs(eps) > TMath::Abs(eps0))
            break;
         eps0 = eps;
      }
      if (s < -TGeoShape::Tolerance())
         continue;
      return TMath::Max(0., s);
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGTorus::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments() + 1;
   nvert = n * (n - 1);
   Bool_t hasrmin = (GetRmin() > 0) ? kTRUE : kFALSE;
   Bool_t hasphi = (GetDphi() < 360) ? kTRUE : kFALSE;
   if (hasrmin)
      nvert *= 2;
   else if (hasphi)
      nvert += 2;
   nsegs = (2 * n - 1) * (n - 1);
   npols = (n - 1) * (n - 1);
   if (hasrmin) {
      nsegs += (2 * n - 1) * (n - 1);
      npols += (n - 1) * (n - 1);
   }
   if (hasphi) {
      nsegs += 2 * (n - 1);
      npols += 2 * (n - 1);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGTorus::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   fBoundingBox.FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments() + 1;
      Int_t nbPnts = n * (n - 1);
      Bool_t hasrmin = (GetRmin() > 0) ? kTRUE : kFALSE;
      Bool_t hasphi = (GetDphi() < 360) ? kTRUE : kFALSE;
      if (hasrmin)
         nbPnts *= 2;
      else if (hasphi)
         nbPnts += 2;

      Int_t nbSegs = (2 * n - 1) * (n - 1);
      Int_t nbPols = (n - 1) * (n - 1);
      if (hasrmin) {
         nbSegs += (2 * n - 1) * (n - 1);
         nbPols += (n - 1) * (n - 1);
      }
      if (hasphi) {
         nbSegs += 2 * (n - 1);
         nbPols += 2 * (n - 1);
      }

      if (buffer.SetRawSizes(nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * nbPols)) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   // TODO: Push down to TGeoShape?? But would have to do raw sizes set first..
   // can rest of TGeoShape be deferred until after
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

