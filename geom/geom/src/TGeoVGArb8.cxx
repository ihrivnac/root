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

// ClassImp(TGeoVGArb8);


////////////////////////////////////////////////////////////////////////////////
/// Default constructor.

TGeoVGArb8::TGeoVGArb8()
  : Base_t("")
{
   fDz = 0;
   for (Int_t i = 0; i < 8; i++) {
      fXY[i][0] = 0.0;
      fXY[i][1] = 0.0;
   }
   SetShapeBit(kGeoArb8);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor. If the array of vertices is not null, this should be
/// in the format : (x0, y0, x1, y1, ... , x7, y7)

TGeoVGArb8::TGeoVGArb8(Double_t dz, Double_t *vertices)
  : Base_t("")
{
   fDz = dz;
   SetShapeBit(kGeoArb8);
   if (vertices) {
      Double_t verticesx[8], verticesy[8];
      for (Int_t i = 0; i < 8; i++) {
         fXY[i][0] = verticesx[i] = vertices[2 * i];
         fXY[i][1] = verticesy[i] = vertices[2 * i + 1];
      }
      // VG shape
      Initialize(verticesx, verticesy, dz);
      ComputeTwist();
      ComputeBBox();
   } else {
      for (Int_t i = 0; i < 8; i++) {
         fXY[i][0] = 0.0;
         fXY[i][1] = 0.0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Named constructor. If the array of vertices is not null, this should be
/// in the format : (x0, y0, x1, y1, ... , x7, y7)

TGeoVGArb8::TGeoVGArb8(const char *name, Double_t dz, Double_t *vertices)
  : Base_t(name)
{
   fDz = dz;
   SetShapeBit(kGeoArb8);
   if (vertices) {
      Double_t verticesx[8], verticesy[8];
      for (Int_t i = 0; i < 8; i++) {
         fXY[i][0] = verticesx[i] = vertices[2 * i];
         fXY[i][1] = verticesy[i] = vertices[2 * i + 1];
      }
      // VG shape
      Initialize(verticesx, verticesy, dz);
      ComputeTwist();
      ComputeBBox();
   } else {
      for (Int_t i = 0; i < 8; i++) {
         fXY[i][0] = 0.0;
         fXY[i][1] = 0.0;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVGArb8::~TGeoVGArb8()
{
   if (fTwist)
      delete[] fTwist;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy twist values from source array

void TGeoVGArb8::CopyTwist(Double_t *twist)
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

void TGeoVGArb8::ComputeBBox()
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

void TGeoVGArb8::ComputeTwist()
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
/// Get twist for segment I in range [0,3]

Double_t TGeoVGArb8::GetTwist(Int_t iseg) const
{
   return (!fTwist || iseg < 0 || iseg > 3) ? 0. : fTwist[iseg];
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this shape along one axis.

TGeoVolume *TGeoVGArb8::Divide(TGeoVolume *voldiv, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Division of an arbitrary trapezoid not implemented");
   return voldiv;
}

////////////////////////////////////////////////////////////////////////////////
/// Get shape range on a given axis.

Double_t TGeoVGArb8::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   if (iaxis == 3) {
      xlo = -fDz;
      xhi = fDz;
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGArb8::GetBoundingCylinder(Double_t *param) const
{
   // first compute rmin/rmax
   Double_t rmaxsq = 0;
   Double_t rsq;
   Int_t i;
   for (i = 0; i < 8; i++) {
      rsq = fXY[i][0] * fXY[i][0] + fXY[i][1] * fXY[i][1];
      rmaxsq = TMath::Max(rsq, rmaxsq);
   }
   param[0] = 0.;     // Rmin
   param[1] = rmaxsq; // Rmax
   param[2] = 0.;     // Phi1
   param[3] = 360.;   // Phi2
}


////////////////////////////////////////////////////////////////////////////////
/// Fills array with n random points located on the surface of indexed facet.
/// The output array must be provided with a length of minimum 3*npoints. Returns
/// true if operation succeeded.
/// Possible index values:
///  - 0 - all facets together
///  - 1 to 6 - facet index from bottom to top Z

Bool_t TGeoVGArb8::GetPointsOnFacet(Int_t /*index*/, Int_t /*npoints*/, Double_t * /* array */) const
{
   return kFALSE;
   /*
      if (index<0 || index>6) return kFALSE;
      if (index==0) {
         // Just generate same number of points on each facet
         Int_t npts = npoints/6.;
         Int_t count = 0;
         for (Int_t ifacet=0; ifacet<6; ifacet++) {
            if (GetPointsOnFacet(ifacet+1, npts, &array[3*count])) count += npts;
            if (ifacet<5) npts = (npoints-count)/(5.-ifacet);
         }
         if (count>0) return kTRUE;
         return kFALSE;
      }
      Double_t z, cf;
      Double_t xmin=TGeoShape::Big();
      Double_t xmax=-xmin;
      Double_t ymin=TGeoShape::Big();
      Double_t ymax=-ymin;
      Double_t dy=0.;
      Double_t poly[8];
      Double_t point[2];
      Int_t i;
      if (index==1 || index==6) {
         z = (index==1)?-fDz:fDz;
         cf = 0.5*(fDz-z)/fDz;
         for (i=0; i<4; i++) {
            poly[2*i]   = fXY[i+4][0]+cf*(fXY[i][0]-fXY[i+4][0]);
            poly[2*i+1] = fXY[i+4][1]+cf*(fXY[i][1]-fXY[i+4][1]);
            xmin = TMath::Min(xmin, poly[2*i]);
            xmax = TMath::Max(xmax, poly[2*i]);
            ymin = TMath::Min(ymin, poly[2*i]);
            ymax = TMath::Max(ymax, poly[2*i]);
         }
      }
      Int_t nshoot = 0;
      Int_t nmiss = 0;
      for (i=0; i<npoints; i++) {
         Double_t *point = &array[3*i];
         switch (surfindex) {
            case 1:
            case 6:
               while (nmiss<1000) {
                  point[0] = xmin + (xmax-xmin)*gRandom->Rndm();
                  point[1] = ymin + (ymax-ymin)*gRandom->Rndm();
               }

      return InsidePolygon(point[0],point[1],poly);
   */
}

////////////////////////////////////////////////////////////////////////////////
/// Prints shape parameters

void TGeoVGArb8::InspectShape() const
{
   printf("*** Shape %s: TGeoVGArb8 ***\n", GetName());
   if (IsTwisted())
      printf("  = TWISTED\n");
   for (Int_t ip = 0; ip < 8; ip++) {
      printf("    point #%i : x=%11.5f y=%11.5f z=%11.5f\n", ip, fXY[ip][0], fXY[ip][1], fDz * ((ip < 4) ? -1 : 1));
   }
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGArb8::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dz       = " << fDz << ";" << std::endl;
   out << "   vert[0]  = " << fXY[0][0] << ";" << std::endl;
   out << "   vert[1]  = " << fXY[0][1] << ";" << std::endl;
   out << "   vert[2]  = " << fXY[1][0] << ";" << std::endl;
   out << "   vert[3]  = " << fXY[1][1] << ";" << std::endl;
   out << "   vert[4]  = " << fXY[2][0] << ";" << std::endl;
   out << "   vert[5]  = " << fXY[2][1] << ";" << std::endl;
   out << "   vert[6]  = " << fXY[3][0] << ";" << std::endl;
   out << "   vert[7]  = " << fXY[3][1] << ";" << std::endl;
   out << "   vert[8]  = " << fXY[4][0] << ";" << std::endl;
   out << "   vert[9]  = " << fXY[4][1] << ";" << std::endl;
   out << "   vert[10] = " << fXY[5][0] << ";" << std::endl;
   out << "   vert[11] = " << fXY[5][1] << ";" << std::endl;
   out << "   vert[12] = " << fXY[6][0] << ";" << std::endl;
   out << "   vert[13] = " << fXY[6][1] << ";" << std::endl;
   out << "   vert[14] = " << fXY[7][0] << ";" << std::endl;
   out << "   vert[15] = " << fXY[7][1] << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGArb8(\"" << GetName() << "\", dz,vert);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Computes intersection points between plane at zpl and non-horizontal edges.

void TGeoVGArb8::SetPlaneVertices(Double_t zpl, Double_t *vertices) const
{
   Double_t cf = 0.5 * (fDz - zpl) / fDz;
   for (Int_t i = 0; i < 4; i++) {
      vertices[2 * i] = fXY[i + 4][0] + cf * (fXY[i][0] - fXY[i + 4][0]);
      vertices[2 * i + 1] = fXY[i + 4][1] + cf * (fXY[i][1] - fXY[i + 4][1]);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Set all arb8 params in one step.
/// param[0] = dz
/// param[1] = x0
/// param[2] = y0
/// ...

void TGeoVGArb8::SetDimensions(Double_t *param)
{
   fDz = param[0];
   Double_t verticesx[8], verticesy[8];
   for (Int_t i = 0; i < 8; i++) {
      fXY[i][0] = verticesx[i] = param[2 * i + 1];
      fXY[i][1] = verticesy[i] = param[2 * i + 2];
   }
   // VG shape
   Initialize(verticesx, verticesy, fDz);
   ComputeTwist();
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates arb8 mesh points

void TGeoVGArb8::SetPoints(Double_t *points) const
{
   for (Int_t i = 0; i < 8; i++) {
      points[3 * i] = fXY[i][0];
      points[3 * i + 1] = fXY[i][1];
      points[3 * i + 2] = (i < 4) ? -fDz : fDz;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Creates arb8 mesh points

void TGeoVGArb8::SetPoints(Float_t *points) const
{
   for (Int_t i = 0; i < 8; i++) {
      points[3 * i] = fXY[i][0];
      points[3 * i + 1] = fXY[i][1];
      points[3 * i + 2] = (i < 4) ? -fDz : fDz;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  Set values for a given vertex.

void TGeoVGArb8::SetVertex(Int_t vnum, Double_t x, Double_t y)
{
   if (vnum < 0 || vnum > 7) {
      Error("SetVertex", "Invalid vertex number");
      return;
   }
   fXY[vnum][0] = x;
   fXY[vnum][1] = y;
   if (vnum == 7) {
      Double_t verticesx[8], verticesy[8];
      for (Int_t i = 0; i < 8; i++) {
         verticesx[i] = fXY[i][0];
         verticesy[i] = fXY[i][1];
      }
      Initialize(verticesx, verticesy, fDz);
      ComputeTwist();
      ComputeBBox();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill size of this 3-D object

void TGeoVGArb8::Sizeof3D() const
{
   TGeoBBox::Sizeof3D();
}

#endif
