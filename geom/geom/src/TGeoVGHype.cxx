// @(#)root/geom:$Id$
// Author: Mihaela Gheata   20/11/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoHype.h"
#include "TGeoVGHype.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGHype);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGHype::TGeoVGHype()
:  Base_t("", 0., 0., 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoHype);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor specifying hyperboloid parameters.

TGeoVGHype::TGeoVGHype(Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
  : Base_t("", rin, rout, stin * TMath::DegToRad(), stout * TMath::DegToRad(), dz)
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetHypeDimensions();
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (dz < 0) {
      SetShapeBit(kGeoRunTimeShape);
   }
   ComputeBBox();
}
////////////////////////////////////////////////////////////////////////////////
/// Constructor specifying parameters and name.

TGeoVGHype::TGeoVGHype(const char *name, Double_t rin, Double_t stin, Double_t rout, Double_t stout, Double_t dz)
  : Base_t(name, rin, rout, stin * TMath::DegToRad(), stout * TMath::DegToRad(), dz)
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetHypeDimensions();
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (dz < 0)
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying a list of parameters
///  - param[0] = dz
///  - param[1] = rin
///  - param[2] = stin
///  - param[3] = rout
///  - param[4] = stout

TGeoVGHype::TGeoVGHype(Double_t *param)
  : Base_t("", param[1], param[3], param[2] * TMath::DegToRad(), param[4] * TMath::DegToRad(), param[0])
{
   SetShapeBit(TGeoShape::kGeoHype);
   SetHypeDimensions();
   // dz<0 can be used to force dz of hyperboloid fit the container volume
   if (param[0] < 0)
      SetShapeBit(kGeoRunTimeShape);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGHype::~TGeoVGHype() {}

////////////////////////////////////////////////////////////////////////////////
/// Compute bounding box of the hyperboloid

void TGeoVGHype::ComputeBBox()
{
   if (GetRmin() < 0.) {
      Warning("ComputeBBox", "Shape %s has invalid rmin=%g ! SET TO 0.", GetName(), GetRmin());
      SetParameters(0., GetRmax(), GetStIn(), GetStOut(), GetDz());
   }
   if ((GetRmin() > GetRmax()) || (GetRmin() * GetRmin() + fTinsq * GetDz() * GetDz() > GetRmax() * GetRmax() + fToutsq * GetDz() * GetDz())) {
      SetShapeBit(kGeoInvalidShape);
      Error("ComputeBBox", "Shape %s hyperbolic surfaces are malformed: rin=%g, stin=%g, rout=%g, stout=%g", GetName(),
            GetRmin(), GetStIn(), GetRmax(), GetStOut());
      return;
   }

   Double_t dx, dy, dz;
   Double_t origin[3];
   dx = dy = TMath::Sqrt(RadiusHypeSq(GetDz(), kFALSE));
   dz = GetDz();
   memset(origin, 0, 3 * sizeof(Double_t));
   fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);
}

////////////////////////////////////////////////////////////////////////////////
/// Cannot divide hyperboloids.

TGeoVolume *TGeoVGHype::Divide(TGeoVolume * /*voldiv*/, const char *divname, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Hyperboloids cannot be divided. Division volume %s not created", divname);
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGHype::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 1: // R
      xlo = GetRmin();
      xhi = TMath::Sqrt(RadiusHypeSq(GetDz(), kFALSE));
      dx = xhi - xlo;
      return dx;
   case 2: // Phi
      xlo = 0;
      xhi = 360;
      dx = 360;
      return dx;
   case 3: // Z
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

void TGeoVGHype::GetBoundingCylinder(Double_t *param) const
{
   param[0] = GetRmin(); // Rmin
   param[0] *= param[0];
   param[1] = TMath::Sqrt(RadiusHypeSq(GetDz(), kFALSE)); // Rmax
   param[1] *= param[1];
   param[2] = 0.;   // Phi1
   param[3] = 360.; // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGHype::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   Double_t dz = GetDz();
   Double_t zmin, zmax;
   if (GetDz() < 0) {
      mother->GetAxisRange(3, zmin, zmax);
      if (zmax < 0)
         return nullptr;
      dz = zmax;
   } else {
      Error("GetMakeRuntimeShape", "Shape %s does not have negative Z range", GetName());
      return nullptr;
   }
   TGeoShape *hype = new TGeoVGHype(GetName(), dz, GetRmax(), GetStOut(), GetRmin(), GetStIn());
   return hype;
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGHype::InspectShape() const
{
   printf("*** Shape %s: TGeoVGHype ***\n", GetName());
   printf("    Rin  = %11.5f\n", GetRmin());
   printf("    sin  = %11.5f\n", GetStIn());
   printf("    Rout = %11.5f\n", GetRmax());
   printf("    sout = %11.5f\n", GetStOut());
   printf("    dz   = %11.5f\n", GetDz());

   printf(" Bounding box:\n");
   fBoundingBox.InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGHype::MakeBuffer3D() const
{
   Int_t n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   Int_t nbPnts = (hasRmin) ? (2 * n * n) : (n * n + 2);
   Int_t nbSegs = (hasRmin) ? (4 * n * n) : (n * (2 * n + 1));
   Int_t nbPols = (hasRmin) ? (2 * n * n) : (n * (n + 1));

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

void TGeoVGHype::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t c = GetBasicColor();
   Int_t i, j, n;
   n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   Int_t irin = 0;
   Int_t irout = (hasRmin) ? (n * n) : 2;
   // Fill segments
   // Case hasRmin:
   //   Inner circles:  [isin = 0], n (per circle) * n ( circles)
   //        iseg = isin+n*i+j , i = 0, n-1   , j = 0, n-1
   //        seg(i=1,n; j=1,n) = [irin+n*i+j] and [irin+n*i+(j+1)%n]
   //   Inner generators: [isgenin = isin+n*n], n (per circle) *(n-1) (slices)
   //        iseg = isgenin + i*n + j, i=0,n-2,  j=0,n-1
   //        seg(i,j) = [irin+n*i+j] and [irin+n*(i+1)+j]
   //   Outer circles:  [isout = isgenin+n*(n-1)], n (per circle) * n ( circles)
   //        iseg = isout + i*n + j , iz = 0, n-1   , j = 0, n-1
   //        seg(i=1,n; j=1,n) = [irout+n*i+j] and [irout+n*i+(j+1)%n]
   //   Outer generators: [isgenout = isout+n*n], n (per circle) *(n-1) (slices)
   //        iseg = isgenout + i*n + j, i=0,n-2,  j=0,n-1
   //        seg(i,j) = [irout+n*i+j] and [irout+n*(i+1)+j]
   //   Lower cap : [islow = isgenout + n*(n-1)], n radial segments
   //        iseg = islow + j,  j=0,n-1
   //        seg(j) = [irin + j] and [irout+j]
   //   Upper cap: [ishi = islow + n], nradial segments
   //        iseg = ishi + j, j=0,n-1
   //        seg[j] = [irin + n*(n-1) + j] and [irout+n*(n-1) + j]
   //
   // Case !hasRmin:
   //   Outer circles: [isout=0], same outer circles (n*n)
   // Outer generators: isgenout = isout + n*n
   //   Lower cap: [islow = isgenout+n*(n-1)], n seg.
   //        iseg = islow + j, j=0,n-1
   //        seg[j] = [irin] and [irout+j]
   //   Upper cap: [ishi = islow +n]
   //        iseg = ishi + j, j=0,n-1
   //        seg[j] = [irin+1] and [irout+n*(n-1) + j]

   Int_t isin = 0;
   Int_t isgenin = (hasRmin) ? (isin + n * n) : 0;
   Int_t isout = (hasRmin) ? (isgenin + n * (n - 1)) : 0;
   Int_t isgenout = isout + n * n;
   Int_t islo = isgenout + n * (n - 1);
   Int_t ishi = islo + n;

   Int_t npt = 0;
   // Fill inner circle segments (n*n)
   if (hasRmin) {
      for (i = 0; i < n; i++) {
         for (j = 0; j < n; j++) {
            npt = 3 * (isin + n * i + j);
            buff.fSegs[npt] = c;
            buff.fSegs[npt + 1] = irin + n * i + j;
            buff.fSegs[npt + 2] = irin + n * i + ((j + 1) % n);
         }
      }
      // Fill inner generators (n*(n-1))
      for (i = 0; i < n - 1; i++) {
         for (j = 0; j < n; j++) {
            npt = 3 * (isgenin + n * i + j);
            buff.fSegs[npt] = c;
            buff.fSegs[npt + 1] = irin + n * i + j;
            buff.fSegs[npt + 2] = irin + n * (i + 1) + j;
         }
      }
   }
   // Fill outer circle segments (n*n)
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         npt = 3 * (isout + n * i + j);
         buff.fSegs[npt] = c;
         buff.fSegs[npt + 1] = irout + n * i + j;
         buff.fSegs[npt + 2] = irout + n * i + ((j + 1) % n);
      }
   }
   // Fill outer generators (n*(n-1))
   for (i = 0; i < n - 1; i++) {
      for (j = 0; j < n; j++) {
         npt = 3 * (isgenout + n * i + j);
         buff.fSegs[npt] = c;
         buff.fSegs[npt + 1] = irout + n * i + j;
         buff.fSegs[npt + 2] = irout + n * (i + 1) + j;
      }
   }
   // Fill lower cap (n)
   for (j = 0; j < n; j++) {
      npt = 3 * (islo + j);
      buff.fSegs[npt] = c;
      buff.fSegs[npt + 1] = irin;
      if (hasRmin)
         buff.fSegs[npt + 1] += j;
      buff.fSegs[npt + 2] = irout + j;
   }
   // Fill upper cap (n)
   for (j = 0; j < n; j++) {
      npt = 3 * (ishi + j);
      buff.fSegs[npt] = c;
      buff.fSegs[npt + 1] = irin + 1;
      if (hasRmin)
         buff.fSegs[npt + 1] += n * (n - 1) + j - 1;
      buff.fSegs[npt + 2] = irout + n * (n - 1) + j;
   }

   // Fill polygons
   // Inner polygons: [ipin = 0] (n-1) slices * n (edges)
   //   ipoly = ipin + n*i + j;  i=0,n-2   j=0,n-1
   //   poly[i,j] = [isin+n*i+j]  [isgenin+i*n+(j+1)%n]  [isin+n*(i+1)+j]  [isgenin+i*n+j]
   // Outer polygons: [ipout = ipin+n*(n-1)]  also (n-1)*n
   //   ipoly = ipout + n*i + j; i=0,n-2   j=0,n-1
   //   poly[i,j] = [isout+n*i+j]  [isgenout+i*n+j]  [isout+n*(i+1)+j]  [isgenout+i*n+(j+1)%n]
   // Lower cap: [iplow = ipout+n*(n-1):  n polygons
   //   ipoly = iplow + j;  j=0,n-1
   //   poly[i=0,j] = [isin+j] [islow+j] [isout+j] [islow+(j+1)%n]
   // Upper cap: [ipup = iplow+n] : n polygons
   //   ipoly = ipup + j;  j=0,n-1
   //   poly[i=n-1, j] = [isin+n*(n-1)+j] [ishi+(j+1)%n] [isout+n*(n-1)+j] [ishi+j]
   //
   // Case !hasRmin:
   // ipin = 0 no inner polygons
   // ipout = 0 same outer polygons
   // Lower cap: iplow = ipout+n*(n-1):  n polygons with 3 segments
   //   poly[i=0,j] = [isout+j] [islow+(j+1)%n] [islow+j]
   // Upper cap: ipup = iplow+n;
   //   poly[i=n-1,j] = [isout+n*(n-1)+j] [ishi+j] [ishi+(j+1)%n]

   Int_t ipin = 0;
   Int_t ipout = (hasRmin) ? (ipin + n * (n - 1)) : 0;
   Int_t iplo = ipout + n * (n - 1);
   Int_t ipup = iplo + n;
   // Inner polygons n*(n-1)
   if (hasRmin) {
      for (i = 0; i < n - 1; i++) {
         for (j = 0; j < n; j++) {
            npt = 6 * (ipin + n * i + j);
            buff.fPols[npt] = c;
            buff.fPols[npt + 1] = 4;
            buff.fPols[npt + 2] = isin + n * i + j;
            buff.fPols[npt + 3] = isgenin + i * n + ((j + 1) % n);
            buff.fPols[npt + 4] = isin + n * (i + 1) + j;
            buff.fPols[npt + 5] = isgenin + i * n + j;
         }
      }
   }
   // Outer polygons n*(n-1)
   for (i = 0; i < n - 1; i++) {
      for (j = 0; j < n; j++) {
         npt = 6 * (ipout + n * i + j);
         buff.fPols[npt] = c;
         buff.fPols[npt + 1] = 4;
         buff.fPols[npt + 2] = isout + n * i + j;
         buff.fPols[npt + 3] = isgenout + i * n + j;
         buff.fPols[npt + 4] = isout + n * (i + 1) + j;
         buff.fPols[npt + 5] = isgenout + i * n + ((j + 1) % n);
      }
   }
   // End caps
   if (hasRmin) {
      for (j = 0; j < n; j++) {
         npt = 6 * (iplo + j);
         buff.fPols[npt] = c + 1;
         buff.fPols[npt + 1] = 4;
         buff.fPols[npt + 2] = isin + j;
         buff.fPols[npt + 3] = islo + j;
         buff.fPols[npt + 4] = isout + j;
         buff.fPols[npt + 5] = islo + ((j + 1) % n);
      }
      for (j = 0; j < n; j++) {
         npt = 6 * (ipup + j);
         buff.fPols[npt] = c + 2;
         buff.fPols[npt + 1] = 4;
         buff.fPols[npt + 2] = isin + n * (n - 1) + j;
         buff.fPols[npt + 3] = ishi + ((j + 1) % n);
         buff.fPols[npt + 4] = isout + n * (n - 1) + j;
         buff.fPols[npt + 5] = ishi + j;
      }
   } else {
      for (j = 0; j < n; j++) {
         npt = 6 * iplo + 5 * j;
         buff.fPols[npt] = c + 1;
         buff.fPols[npt + 1] = 3;
         buff.fPols[npt + 2] = isout + j;
         buff.fPols[npt + 3] = islo + ((j + 1) % n);
         buff.fPols[npt + 4] = islo + j;
      }
      for (j = 0; j < n; j++) {
         npt = 6 * iplo + 5 * (n + j);
         buff.fPols[npt] = c + 2;
         buff.fPols[npt + 1] = 3;
         buff.fPols[npt + 2] = isout + n * (n - 1) + j;
         buff.fPols[npt + 3] = ishi + j;
         buff.fPols[npt + 4] = ishi + ((j + 1) % n);
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return StIn

Double_t TGeoVGHype::GetStIn() const
{
   return Base_t::GetStIn() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Return StOut

Double_t TGeoVGHype::GetStOut() const
{
   return Base_t::GetStOut() * TMath::RadToDeg();
}

////////////////////////////////////////////////////////////////////////////////
/// Compute r^2 = x^2 + y^2 at a given z coordinate, for either inner or outer hyperbolas.

Double_t TGeoVGHype::RadiusHypeSq(Double_t z, Bool_t inner) const
{
   Double_t r0, tsq;
   if (inner) {
      r0 = GetRmin();
      tsq = fTinsq;
   } else {
      r0 = GetRmax();
      tsq = fToutsq;
   }
   return (r0 * r0 + tsq * z * z);
}

////////////////////////////////////////////////////////////////////////////////
/// Compute z^2 at a given  r^2, for either inner or outer hyperbolas.

Double_t TGeoVGHype::ZHypeSq(Double_t r, Bool_t inner) const
{
   Double_t r0, tsq;
   if (inner) {
      r0 = GetRmin();
      tsq = fTinsq;
   } else {
      r0 = GetRmax();
      tsq = fToutsq;
   }
   if (TMath::Abs(tsq) < TGeoShape::Tolerance())
      return TGeoShape::Big();
   return ((r * r - r0 * r0) / tsq);
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGHype::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   rin   = " << GetRmin() << ";" << std::endl;
   out << "   stin  = " << GetStIn() << ";" << std::endl;
   out << "   rout  = " << GetRmax() << ";" << std::endl;
   out << "   stout = " << GetStOut() << ";" << std::endl;
   out << "   dz    = " << GetDz() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGHype(\"" << GetName() << "\",rin,stin,rout,stout,dz);"
       << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the hyperboloid.

void TGeoVGHype::SetHypeDimensions()
{
   fTin = TMath::Tan(Base_t::GetStIn());
   fTinsq = fTin * fTin;
   fTout = TMath::Tan(Base_t::GetStOut());
   fToutsq = fTout * fTout;
   if ((GetRmin() == 0) && (GetStIn() == 0))
      SetShapeBit(kGeoRSeg, kTRUE);
   else
      SetShapeBit(kGeoRSeg, kFALSE);
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoVGHype::SetPoints(Double_t *points) const
{
   Double_t z, dz, r;
   Int_t i, j, n;
   if (!points)
      return;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   dz = 2. * GetDz() / (n - 1);

   Int_t indx = 0;

   if (HasInner()) {
      // Inner surface points
      for (i = 0; i < n; i++) {
         z = -GetDz() + i * dz;
         r = TMath::Sqrt(RadiusHypeSq(z, kTRUE));
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx++] = r * TMath::Cos(phi);
            points[indx++] = r * TMath::Sin(phi);
            points[indx++] = z;
         }
      }
   } else {
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -GetDz();
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = GetDz();
   }
   // Outer surface points
   for (i = 0; i < n; i++) {
      z = -GetDz() + i * dz;
      r = TMath::Sqrt(RadiusHypeSq(z, kFALSE));
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = r * TMath::Cos(phi);
         points[indx++] = r * TMath::Sin(phi);
         points[indx++] = z;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create tube mesh points

void TGeoVGHype::SetPoints(Float_t *points) const
{
   Double_t z, dz, r;
   Int_t i, j, n;
   if (!points)
      return;
   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   dz = 2. * GetDz() / (n - 1);

   Int_t indx = 0;

   if (HasInner()) {
      // Inner surface points
      for (i = 0; i < n; i++) {
         z = -GetDz() + i * dz;
         r = TMath::Sqrt(RadiusHypeSq(z, kTRUE));
         for (j = 0; j < n; j++) {
            phi = j * dphi * TMath::DegToRad();
            points[indx++] = r * TMath::Cos(phi);
            points[indx++] = r * TMath::Sin(phi);
            points[indx++] = z;
         }
      }
   } else {
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = -GetDz();
      points[indx++] = 0.;
      points[indx++] = 0.;
      points[indx++] = GetDz();
   }
   // Outer surface points
   for (i = 0; i < n; i++) {
      z = -GetDz() + i * dz;
      r = TMath::Sqrt(RadiusHypeSq(z, kFALSE));
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         points[indx++] = r * TMath::Cos(phi);
         points[indx++] = r * TMath::Sin(phi);
         points[indx++] = z;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGHype::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t n = gGeoManager->GetNsegments();
   Bool_t hasRmin = HasInner();
   nvert = (hasRmin) ? (2 * n * n) : (n * n + 2);
   nsegs = (hasRmin) ? (4 * n * n) : (n * (2 * n + 1));
   npols = (hasRmin) ? (2 * n * n) : (n * (n + 1));
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGHype::GetNmeshVertices() const
{
   Int_t n = gGeoManager->GetNsegments();
   Int_t numPoints = (HasRmin()) ? (2 * n * n) : (n * n + 2);
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGHype::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGHype::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   fBoundingBox.FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t n = gGeoManager->GetNsegments();
      Bool_t hasRmin = HasInner();
      Int_t nbPnts = (hasRmin) ? (2 * n * n) : (n * n + 2);
      Int_t nbSegs = (hasRmin) ? (4 * n * n) : (n * (2 * n + 1));
      Int_t nbPols = (hasRmin) ? (2 * n * n) : (n * (n + 1));
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
