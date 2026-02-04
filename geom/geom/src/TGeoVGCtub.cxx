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

#include "TGeoTube.h"
#include "TGeoVGTube.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGCtub);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGCtub::TGeoVGCtub()
   : Base_t("", 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
{
   // default ctor
   SetShapeBit(kGeoCtub);
   fNlow[0] = fNlow[1] = fNhigh[0] = fNhigh[1] = 0.;
   fNlow[2] = -1;
   fNhigh[2] = 1;
}

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

TGeoVGCtub::TGeoVGCtub(Double_t *params)
  : Base_t("", 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
{
   SetCtubDimensions(params[0], params[1], params[2], params[3], params[4], params[5], params[6], params[7], params[8],
                     params[9], params[10]);
   SetShapeBit(kGeoCtub);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGCtub::~TGeoVGCtub() {}

// ////////////////////////////////////////////////////////////////////////////////
// /// Init frequently used trigonometric values

void TGeoVGCtub::InitTrigonometry()
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
/// compute bounding box of the tube segment
/// (copied) from TGeoVGTubeSeg::ComputeBBox()

void TGeoVGCtub::ComputeTubeSegBBox()
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
/// compute minimum bounding box of the ctub

void TGeoVGCtub::ComputeBBox()
{
   ComputeTubeSegBBox();

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
      Double_t dphi = Base_t::dphi() * TMath::RadToDeg();
      if (dphi < 0)
         dphi += 360.;
      Double_t ddp = phi_low - sphi() * TMath::RadToDeg();
      if (ddp < 0)
         ddp += 360.;
      if (ddp <= dphi) {
         xc = rmin() * TMath::Cos(phi_low * TMath::DegToRad());
         yc = rmin() * TMath::Sin(phi_low * TMath::DegToRad());
         z1 = GetZcoord(xc, yc, -Base_t::z());
         xc = rmax() * TMath::Cos(phi_low * TMath::DegToRad());
         yc = rmax() * TMath::Sin(phi_low * TMath::DegToRad());
         z1 = TMath::Min(z1, GetZcoord(xc, yc, -Base_t::z()));
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
      Double_t dphi = Base_t::dphi() * TMath::RadToDeg();;
      if (dphi < 0)
         dphi += 360.;
      Double_t ddp = phi_hi - sphi() * TMath::RadToDeg();
      if (ddp < 0)
         ddp += 360.;
      if (ddp <= dphi) {
         xc = rmin() * TMath::Cos(phi_hi * TMath::DegToRad());
         yc = rmin() * TMath::Sin(phi_hi * TMath::DegToRad());
         z1 = GetZcoord(xc, yc, Base_t::z());
         xc = rmax() * TMath::Cos(phi_hi * TMath::DegToRad());
         yc = rmax() * TMath::Sin(phi_hi * TMath::DegToRad());
         z1 = TMath::Max(z1, GetZcoord(xc, yc, Base_t::z()));
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
   z[0] = GetZcoord(xc, yc, -Base_t::z());
   z[4] = GetZcoord(xc, yc, Base_t::z());

   xc = rmin() * fC2;
   yc = rmin() * fS2;
   z[1] = GetZcoord(xc, yc, -Base_t::z());
   z[5] = GetZcoord(xc, yc, Base_t::z());

   xc = rmax() * fC1;
   yc = rmax() * fS1;
   z[2] = GetZcoord(xc, yc, -Base_t::z());
   z[6] = GetZcoord(xc, yc, Base_t::z());

   xc = rmax() * fC2;
   yc = rmax() * fS2;
   z[3] = GetZcoord(xc, yc, -Base_t::z());
   z[7] = GetZcoord(xc, yc, Base_t::z());

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
      xlo = sphi() * TMath::RadToDeg();
      xhi = (sphi() + dphi()) * TMath::RadToDeg();
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
      dz = ((TGeoVGCtub *)mother)->GetDz();
   if (Base_t::rmin() < 0)
      rmin = ((TGeoVGCtub *)mother)->GetRmin();
   if ((Base_t::rmax() < 0) || (Base_t::rmax() <= Base_t::rmin()))
      rmax = ((TGeoVGCtub *)mother)->GetRmax();

   return (new TGeoVGCtub(rmin, rmax, dz, sphi() * TMath::RadToDeg(), (sphi() + dphi()) * TMath::RadToDeg(), fNlow[0], fNlow[1], fNlow[2], fNhigh[0], fNhigh[1], fNhigh[2]));
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

void TGeoVGCtub::SetTubsDimensions(Double_t rmin, Double_t rmax, Double_t dz, Double_t phiStart, Double_t phiEnd)
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
   phi1 = sphi() * TMath::RadToDeg();
   phi2 = (sphi() + Base_t::dphi()) * TMath::RadToDeg();
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
   phi1 = sphi() * TMath::RadToDeg();
   phi2 = (sphi() + Base_t::dphi()) * TMath::RadToDeg();
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
      buffer.fPhiMin = sphi() * TMath::RadToDeg();
      buffer.fPhiMax = (sphi() + dphi()) * TMath::RadToDeg();

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

#endif // ROOT_USE_VECGEOM_SOLIDS
