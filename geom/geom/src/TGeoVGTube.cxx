// @(#)root/geom:$Id$

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
#include "TGeoTube.h"   // For tubeseg, ctubs
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <VecGeom/base/Math.h>

#include <iostream>

// ClassImp(TGeoVGTube);

// ////////////////////////////////////////////////////////////////////////////////
// /// Default constructor

TGeoVGTube::TGeoVGTube()
 : Base_t("", 0., 0., 0., 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoTube);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGTube::TGeoVGTube(Double_t rmin, Double_t rmax, Double_t dz)
 : Base_t("", rmin, rmax, dz, 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoTube);
   // SetTubeDimensions(rmin, rmax, dz);
   if ((dz < 0) || (rmin < 0) || (rmax < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   ComputeBBox();
}
////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius

TGeoVGTube::TGeoVGTube(const char *name, Double_t rmin, Double_t rmax, Double_t dz)
  : Base_t(name, rmin, rmax, dz, 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoTube);
   if ((dz < 0) || (rmin < 0) || (rmax < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
///  - param[0] = Rmin
///  - param[1] = Rmax
///  - param[2] = dz

TGeoVGTube::TGeoVGTube(Double_t *param)
  : Base_t("", param[0], param[1], param[2], 0., vecgeom::kTwoPi)
{
   SetShapeBit(TGeoShape::kGeoTube);
   if ((param[0] < 0) || (param[1] < 0) || (param[2] < 0))
      SetShapeBit(kGeoRunTimeShape);
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the tube

void TGeoVGTube::ComputeBBox()
{
// get parameters from vecgeom

   fDX = fDY = rmax();
   fDZ = z();
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Computes capacity of the shape in [length^3] 

Double_t TGeoVGTube::Capacity(Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t capacity = 2. * TMath::Pi() * (rmax * rmax - rmin * rmin) * dz;
   return capacity;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute normal to closest surface from POINT.

void TGeoVGTube::ComputeNormalS(const Double_t *point, const Double_t *dir, Double_t *norm, Double_t /*rmin*/,
                              Double_t /*rmax*/, Double_t /*dz*/)
{
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
/// Static: Compute distance from inside point to surface of the tube (static)
/// Boundary safe algorithm.
/// compute distance to surface
/// Do Z

Double_t
TGeoVGTube::DistFromInsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t sz = TGeoShape::Big();
   if (dir[2]) {
      sz = (TMath::Sign(dz, dir[2]) - point[2]) / dir[2];
      if (sz <= 0)
         return 0.0;
   }
   // Do R
   Double_t nsq = dir[0] * dir[0] + dir[1] * dir[1];
   if (TMath::Abs(nsq) < TGeoShape::Tolerance())
      return sz;
   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
   Double_t b, d;
   // inner cylinder
   if (rmin > 0) {
      // Protection in case point is actually outside the tube
      if (rsq <= rmin * rmin + TGeoShape::Tolerance()) {
         if (rdotn < 0)
            return 0.0;
      } else {
         if (rdotn < 0) {
            DistToTube(rsq, nsq, rdotn, rmin, b, d);
            if (d > 0) {
               Double_t sr = -b - d;
               if (sr > 0)
                  return TMath::Min(sz, sr);
            }
         }
      }
   }
   // outer cylinder
   if (rsq >= rmax * rmax - TGeoShape::Tolerance()) {
      if (rdotn >= 0)
         return 0.0;
   }
   DistToTube(rsq, nsq, rdotn, rmax, b, d);
   if (d > 0) {
      Double_t sr = -b + d;
      if (sr > 0)
         return TMath::Min(sz, sr);
   }
   return 0.;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to compute distance from outside point to a tube with given parameters
/// Boundary safe algorithm.
/// check Z planes

Double_t
TGeoVGTube::DistFromOutsideS(const Double_t *point, const Double_t *dir, Double_t rmin, Double_t rmax, Double_t dz)
{
   Double_t xi, yi, zi;
   Double_t rmaxsq = rmax * rmax;
   Double_t rminsq = rmin * rmin;
   zi = dz - TMath::Abs(point[2]);
   Bool_t in = kFALSE;
   Bool_t inz = (zi < 0) ? kFALSE : kTRUE;
   if (!inz) {
      if (point[2] * dir[2] >= 0)
         return TGeoShape::Big();
      Double_t s = -zi / TMath::Abs(dir[2]);
      xi = point[0] + s * dir[0];
      yi = point[1] + s * dir[1];
      Double_t r2 = xi * xi + yi * yi;
      if ((rminsq <= r2) && (r2 <= rmaxsq))
         return s;
   }

   Double_t rsq = point[0] * point[0] + point[1] * point[1];
   // check outer cyl. surface
   Double_t nsq = dir[0] * dir[0] + dir[1] * dir[1];
   Double_t rdotn = point[0] * dir[0] + point[1] * dir[1];
   Double_t b, d;
   Bool_t inrmax = kFALSE;
   Bool_t inrmin = kFALSE;
   if (rsq <= rmaxsq + TGeoShape::Tolerance())
      inrmax = kTRUE;
   if (rsq >= rminsq - TGeoShape::Tolerance())
      inrmin = kTRUE;
   in = inz & inrmin & inrmax;
   // If inside, we are most likely on a boundary within machine precision.
   if (in) {
      Bool_t checkout = kFALSE;
      Double_t r = TMath::Sqrt(rsq);
      if (zi < rmax - r) {
         if ((TGeoShape::IsSameWithinTolerance(rmin, 0)) || (zi < r - rmin)) {
            if (point[2] * dir[2] < 0)
               return 0.0;
            return TGeoShape::Big();
         }
      }
      if ((rmaxsq - rsq) < (rsq - rminsq))
         checkout = kTRUE;
      if (checkout) {
         if (rdotn >= 0)
            return TGeoShape::Big();
         return 0.0;
      }
      if (TGeoShape::IsSameWithinTolerance(rmin, 0))
         return 0.0;
      if (rdotn >= 0)
         return 0.0;
      // Ray exiting rmin -> check (+) solution for inner tube
      if (TMath::Abs(nsq) < TGeoShape::Tolerance())
         return TGeoShape::Big();
      DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d > 0) {
         Double_t s = -b + d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz)
               return s;
         }
      }
      return TGeoShape::Big();
   }
   // Check outer cylinder (only r>rmax has to be considered)
   if (TMath::Abs(nsq) < TGeoShape::Tolerance())
      return TGeoShape::Big();
   if (!inrmax) {
      DistToTube(rsq, nsq, rdotn, rmax, b, d);
      if (d > 0) {
         Double_t s = -b - d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz)
               return s;
         }
      }
   }
   // check inner cylinder
   if (rmin > 0) {
      DistToTube(rsq, nsq, rdotn, rmin, b, d);
      if (d > 0) {
         Double_t s = -b + d;
         if (s > 0) {
            zi = point[2] + s * dir[2];
            if (TMath::Abs(zi) <= dz)
               return s;
         }
      }
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// Static method computing the distance to a tube with given radius, starting from
/// POINT along DIR director cosines. The distance is computed as :
///    RSQ   = point[0]*point[0]+point[1]*point[1]
///    NSQ   = dir[0]*dir[0]+dir[1]*dir[1]  ---> should NOT be 0 !!!
///    RDOTN = point[0]*dir[0]+point[1]*dir[1]
/// The distance can be computed as :
///    D = -B +/- DELTA
/// where DELTA.GT.0 and D.GT.0

void TGeoVGTube::DistToTube(Double_t rsq, Double_t nsq, Double_t rdotn, Double_t radius, Double_t &b, Double_t &delta)
{
   Double_t t1 = 1. / nsq;
   Double_t t3 = rsq - (radius * radius);
   b = t1 * rdotn;
   Double_t c = t1 * t3;
   delta = b * b - c;
   if (delta > 0) {
      delta = TMath::Sqrt(delta);
   } else {
      delta = -1;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this tube shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume in case of Z divisions. For radial division
/// creates all volumes with different shapes and returns pointer to volume that
/// was divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *
TGeoVGTube::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
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
      finder = new TGeoPatternCylR(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         shape = new TGeoVGTube(start + id * step, start + (id + 1) * step, z());
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "R";
         voldiv->AddNodeOffset(vol, id, 0, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   case 2: //---                Phi division
      finder = new TGeoPatternCylPhi(voldiv, ndiv, start, end);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      // TODO Change to VG
      shape = new TGeoTubeSeg(rmin(), rmax(), z(), -step / 2, step / 2);
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
      finder = new TGeoPatternZ(voldiv, ndiv, start, start + ndiv * step);
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      shape = new TGeoVGTube(rmin(), rmax(), step / 2);
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
/// Returns name of axis IAXIS.

const char *TGeoVGTube::GetAxisName(Int_t iaxis) const
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

Double_t TGeoVGTube::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
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
      xlo = 0;
      xhi = 360;
      dx = 360;
      return dx;
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
/// is the following : Rmin, Rmax, Phi1, Phi2, dZ

void TGeoVGTube::GetBoundingCylinder(Double_t *param) const
{
   param[0] = rmin(); // Rmin
   param[0] *= param[0];
   param[1] = rmax(); // Rmax
   param[1] *= param[1];
   param[2] = 0.;   // Phi1
   param[3] = 360.; // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGTube::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   Double_t rmin, rmax, dz;
   Double_t xmin, xmax;
   rmin = Base_t::rmin();
   rmax = Base_t::rmax();
   dz = z();
   if (z() < 0) {
      mother->GetAxisRange(3, xmin, xmax);
      if (xmax < 0)
         return nullptr;
      dz = xmax;
   }
   mother->GetAxisRange(1, xmin, xmax);
   if (Base_t::rmin() < 0) {
      if (xmin < 0)
         return nullptr;
      rmin = xmin;
   }
   if (Base_t::rmax() < 0) {
      if (xmax <= 0)
         return nullptr;
      rmax = xmax;
   }

   return (new TGeoVGTube(GetName(), rmin, rmax, dz));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGTube::InspectShape() const
{
   printf("*** Shape %s: TGeoVGTube ***\n", GetName());
   printf("    Rmin = %11.5f\n", rmin());
   printf("    Rmax = %11.5f\n", rmax());
   printf("    dz   = %11.5f\n", z());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGTube::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t nbPnts = 4 * n;
   Int_t nbSegs = 8 * n;
   Int_t nbPols = 4 * n;
   if (!HasRmin()) {
      nbPnts = 2 * (n + 1);
      nbSegs = 5 * n;
      nbPols = 3 * n;
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
/// Fill TBuffer3D structure for segments and polygons.

void TGeoVGTube::SetSegsAndPols(TBuffer3D &buffer) const
{
   Int_t i, j, indx;
   Int_t n = gGeoManager->GetNsegments();
   Int_t c = (((buffer.fColor) % 8) - 1) * 4;
   if (c < 0)
      c = 0;

   if (HasRmin()) {
      // circle segments:
      // lower rmin circle: i=0, (0, n-1)
      // lower rmax circle: i=1, (n, 2n-1)
      // upper rmin circle: i=2, (2n, 3n-1)
      // upper rmax circle: i=1, (3n, 4n-1)
      for (i = 0; i < 4; i++) {
         for (j = 0; j < n; j++) {
            indx = 3 * (i * n + j);
            buffer.fSegs[indx] = c;
            buffer.fSegs[indx + 1] = i * n + j;
            buffer.fSegs[indx + 2] = i * n + (j + 1) % n;
         }
      }
      // Z-parallel segments
      // inner: i=4, (4n, 5n-1)
      // outer: i=5, (5n, 6n-1)
      for (i = 4; i < 6; i++) {
         for (j = 0; j < n; j++) {
            indx = 3 * (i * n + j);
            buffer.fSegs[indx] = c + 1;
            buffer.fSegs[indx + 1] = (i - 4) * n + j;
            buffer.fSegs[indx + 2] = (i - 2) * n + j;
         }
      }
      // Radial segments
      // lower: i=6, (6n, 7n-1)
      // upper: i=7, (7n, 8n-1)
      for (i = 6; i < 8; i++) {
         for (j = 0; j < n; j++) {
            indx = 3 * (i * n + j);
            buffer.fSegs[indx] = c;
            buffer.fSegs[indx + 1] = 2 * (i - 6) * n + j;
            buffer.fSegs[indx + 2] = (2 * (i - 6) + 1) * n + j;
         }
      }
      // Polygons
      i = 0;
      // Inner lateral (0, n-1)
      for (j = 0; j < n; j++) {
         indx = 6 * (i * n + j);
         buffer.fPols[indx] = c;
         buffer.fPols[indx + 1] = 4;
         buffer.fPols[indx + 2] = j;
         buffer.fPols[indx + 3] = 4 * n + (j + 1) % n;
         buffer.fPols[indx + 4] = 2 * n + j;
         buffer.fPols[indx + 5] = 4 * n + j;
      }
      i = 1;
      // Outer lateral (n,2n-1)
      for (j = 0; j < n; j++) {
         indx = 6 * (i * n + j);
         buffer.fPols[indx] = c + 1;
         buffer.fPols[indx + 1] = 4;
         buffer.fPols[indx + 2] = n + j;
         buffer.fPols[indx + 3] = 5 * n + j;
         buffer.fPols[indx + 4] = 3 * n + j;
         buffer.fPols[indx + 5] = 5 * n + (j + 1) % n;
      }
      i = 2;
      // lower disc (2n, 3n-1)
      for (j = 0; j < n; j++) {
         indx = 6 * (i * n + j);
         buffer.fPols[indx] = c;
         buffer.fPols[indx + 1] = 4;
         buffer.fPols[indx + 2] = j;
         buffer.fPols[indx + 3] = 6 * n + j;
         buffer.fPols[indx + 4] = n + j;
         buffer.fPols[indx + 5] = 6 * n + (j + 1) % n;
      }
      i = 3;
      // upper disc (3n, 4n-1)
      for (j = 0; j < n; j++) {
         indx = 6 * (i * n + j);
         buffer.fPols[indx] = c;
         buffer.fPols[indx + 1] = 4;
         buffer.fPols[indx + 2] = 2 * n + j;
         buffer.fPols[indx + 3] = 7 * n + (j + 1) % n;
         buffer.fPols[indx + 4] = 3 * n + j;
         buffer.fPols[indx + 5] = 7 * n + j;
      }
      return;
   }
   // Rmin=0 tubes
   // circle segments
   // lower rmax circle: i=0, (0, n-1)
   // upper rmax circle: i=1, (n, 2n-1)
   for (i = 0; i < 2; i++) {
      for (j = 0; j < n; j++) {
         indx = 3 * (i * n + j);
         buffer.fSegs[indx] = c;
         buffer.fSegs[indx + 1] = 2 + i * n + j;
         buffer.fSegs[indx + 2] = 2 + i * n + (j + 1) % n;
      }
   }
   // Z-parallel segments (2n,3n-1)
   for (j = 0; j < n; j++) {
      indx = 3 * (2 * n + j);
      buffer.fSegs[indx] = c + 1;
      buffer.fSegs[indx + 1] = 2 + j;
      buffer.fSegs[indx + 2] = 2 + n + j;
   }
   // Radial segments
   // Lower circle: i=3, (3n,4n-1)
   // Upper circle: i=4, (4n,5n-1)
   for (i = 3; i < 5; i++) {
      for (j = 0; j < n; j++) {
         indx = 3 * (i * n + j);
         buffer.fSegs[indx] = c;
         buffer.fSegs[indx + 1] = i - 3;
         buffer.fSegs[indx + 2] = 2 + (i - 3) * n + j;
      }
   }
   // Polygons
   // lateral (0,n-1)
   for (j = 0; j < n; j++) {
      indx = 6 * j;
      buffer.fPols[indx] = c + 1;
      buffer.fPols[indx + 1] = 4;
      buffer.fPols[indx + 2] = j;
      buffer.fPols[indx + 3] = 2 * n + j;
      buffer.fPols[indx + 4] = n + j;
      buffer.fPols[indx + 5] = 2 * n + (j + 1) % n;
   }
   // bottom triangles (n,2n-1)
   for (j = 0; j < n; j++) {
      indx = 6 * n + 5 * j;
      buffer.fPols[indx] = c;
      buffer.fPols[indx + 1] = 3;
      buffer.fPols[indx + 2] = j;
      buffer.fPols[indx + 3] = 3 * n + (j + 1) % n;
      buffer.fPols[indx + 4] = 3 * n + j;
   }
   // top triangles (2n,3n-1)
   for (j = 0; j < n; j++) {
      indx = 6 * n + 5 * n + 5 * j;
      buffer.fPols[indx] = c;
      buffer.fPols[indx + 1] = 3;
      buffer.fPols[indx + 2] = n + j;
      buffer.fPols[indx + 3] = 4 * n + j;
      buffer.fPols[indx + 4] = 4 * n + (j + 1) % n;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Static: computes the closest distance from given point to this shape, according
/// to option. The matching point on the shape is stored in spoint.

Double_t TGeoVGTube::SafetyS(const Double_t *point, Bool_t in, Double_t rmin, Double_t rmax, Double_t dz, Int_t skipz)
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
   saf[1] = (rmin > 1E-10) ? (r - rmin) : TGeoShape::Big();
   saf[2] = rmax - r;
   //   printf("saf0=%g saf1=%g saf2=%g in=%d skipz=%d\n", saf[0],saf[1],saf[2],in,skipz);
   if (in)
      return saf[TMath::LocMin(3, saf)];
   for (Int_t i = 0; i < 3; i++)
      saf[i] = -saf[i];
   return saf[TMath::LocMax(3, saf)];
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGTube::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rmin = " << rmin() << ";" << std::endl;
   out << "   rmax = " << rmax() << ";" << std::endl;
   out << "   dz   = " << z() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGTube(\"" << GetName() << "\",rmin,rmax,dz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tube dimensions.

void TGeoVGTube::SetTubeDimensions(Double_t rmin, Double_t rmax, Double_t dz)
{
   SetRMin(rmin);
   SetRMin(rmax);
   SetDz(dz);
   if (rmin > 0 && rmax > 0 && rmin >= rmax)
      Error("SetTubeDimensions", "In shape %s wrong rmin=%g rmax=%g", GetName(), rmin, rmax);
}

////////////////////////////////////////////////////////////////////////////////
/// Set tube dimensions starting from a list.

void TGeoVGTube::SetDimensions(Double_t *param)
{
   Double_t rmin = param[0];
   Double_t rmax = param[1];
   Double_t dz = param[2];
   SetTubeDimensions(rmin, rmax, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the line segments of the shape mesh.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation is implemented.

Bool_t TGeoVGTube::GetPointsOnSegments(Int_t npoints, Double_t *array) const
{
   if (npoints > (npoints / 2) * 2) {
      Error("GetPointsOnSegments", "Npoints must be even number");
      return kFALSE;
   }
   Int_t nc = 0;
   if (HasRmin())
      nc = (Int_t)TMath::Sqrt(0.5 * npoints);
   else
      nc = (Int_t)TMath::Sqrt(1. * npoints);
   Double_t dphi = TMath::TwoPi() / nc;
   Double_t phi = 0;
   Int_t ntop = 0;
   if (HasRmin())
      ntop = npoints / 2 - nc * (nc - 1);
   else
      ntop = npoints - nc * (nc - 1);
   Double_t dz = 2 * z() / (nc - 1);
   Double_t z = 0;
   Int_t icrt = 0;
   Int_t nphi = nc;
   // loop z sections
   for (Int_t i = 0; i < nc; i++) {
      if (i == (nc - 1))
         nphi = ntop;
      z = -Base_t::z() + i * dz;
      // loop points on circle sections
      for (Int_t j = 0; j < nphi; j++) {
         phi = j * dphi;
         if (HasRmin()) {
            array[icrt++] = rmin() * TMath::Cos(phi);
            array[icrt++] = rmin() * TMath::Sin(phi);
            array[icrt++] = z;
         }
         array[icrt++] = rmax() * TMath::Cos(phi);
         array[icrt++] = rmax() * TMath::Sin(phi);
         array[icrt++] = z;
      }
   }
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoVGTube::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   dz = z();
   Int_t indx = 0;
   if (points) {
      if (HasRmin()) {
         // 4*n points
         // (0,n-1) lower rmin circle
         // (2n, 3n-1) upper rmin circle
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
            indx++;
            points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
            indx++;
            points[indx + 6 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
         // (n, 2n-1) lower rmax circle
         // (3n, 4n-1) upper rmax circle
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
            indx++;
            points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
            indx++;
            points[indx + 6 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
      } else {
         // centers of lower/upper circles (0,1)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -dz;
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = dz;
         // lower rmax circle (2, 2+n-1)
         // upper rmax circle (2+n, 2+2n-1)
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 3 * n] = points[indx] = rmax() * TMath::Cos(phi);
            indx++;
            points[indx + 3 * n] = points[indx] = rmax() * TMath::Sin(phi);
            indx++;
            points[indx + 3 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoVGTube::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   dz = z();
   Int_t indx = 0;
   if (points) {
      if (HasRmin()) {
         // 4*n points
         // (0,n-1) lower rmin circle
         // (2n, 3n-1) upper rmin circle
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 6 * n] = points[indx] = rmin() * TMath::Cos(phi);
            indx++;
            points[indx + 6 * n] = points[indx] = rmin() * TMath::Sin(phi);
            indx++;
            points[indx + 6 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
         // (n, 2n-1) lower rmax circle
         // (3n, 4n-1) upper rmax circle
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 6 * n] = points[indx] = rmax() * TMath::Cos(phi);
            indx++;
            points[indx + 6 * n] = points[indx] = rmax() * TMath::Sin(phi);
            indx++;
            points[indx + 6 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
      } else {
         // centers of lower/upper circles (0,1)
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = -dz;
         points[indx++] = 0.;
         points[indx++] = 0.;
         points[indx++] = dz;
         // lower rmax circle (2, 2+n-1)
         // upper rmax circle (2+n, 2+2n-1)
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx + 3 * n] = points[indx] = rmax() * TMath::Cos(phi);
            indx++;
            points[indx + 3 * n] = points[indx] = rmax() * TMath::Sin(phi);
            indx++;
            points[indx + 3 * n] = dz;
            points[indx] = -dz;
            indx++;
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGTube::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = n * 4;
   if (!HasRmin())
      numPoints = 2 * (n + 1);
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGTube::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments();
   nvert = n * 4;
   nsegs = n * 8;
   npols = n * 4;
   if (!HasRmin()) {
      nvert = 2 * (n + 1);
      nsegs = 5 * n;
      npols = 3 * n;
   } else {
      nvert = n * 4;
      nsegs = n * 8;
      npols = n * 4;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGTube::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGTube::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3DTube buffer;
   TGeoBBox::FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kShapeSpecific) {
      buffer.fRadiusInner = rmin();
      buffer.fRadiusOuter = rmax();
      buffer.fHalfLength = z();
      buffer.SetSectionsValid(TBuffer3D::kShapeSpecific);
   }
   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments();
      Int_t nbPnts = 4 * n;
      Int_t nbSegs = 8 * n;
      Int_t nbPols = 4 * n;
      if (!HasRmin()) {
         nbPnts = 2 * (n + 1);
         nbSegs = 5 * n;
         nbPols = 3 * n;
      }
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

#endif // ROOT_USE_VECGEOM_SOLIDS
