// @(#)root/geom:$Id$
// Author: Mihaela Gheata   24/01/04

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#include "TGeoVGXtru.h"
#include "TGeoVGXtru.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>

#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include "TVirtualGeoPainter.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoPolygon.h"

// ClassImp(TGeoVGXtru);

////////////////////////////////////////////////////////////////////////////////
/// Constructor.

TGeoVGXtru::ThreadData_t::ThreadData_t() : fSeg(0), fIz(0), fXc(nullptr), fYc(nullptr), fPoly(nullptr) {}

////////////////////////////////////////////////////////////////////////////////
/// Destructor.

TGeoVGXtru::ThreadData_t::~ThreadData_t()
{
   delete[] fXc;
   delete[] fYc;
   delete fPoly;
}

////////////////////////////////////////////////////////////////////////////////

TGeoVGXtru::ThreadData_t &TGeoVGXtru::GetThreadData() const
{
   if (!fThreadSize)
      ((TGeoVGXtru *)this)->CreateThreadData(1);
   Int_t tid = TGeoManager::ThreadId();
   return *fThreadData[tid];
}

////////////////////////////////////////////////////////////////////////////////

void TGeoVGXtru::ClearThreadData() const
{
   std::lock_guard<std::mutex> guard(fMutex);
   std::vector<ThreadData_t *>::iterator i = fThreadData.begin();
   while (i != fThreadData.end()) {
      delete *i;
      ++i;
   }
   fThreadData.clear();
   fThreadSize = 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Create thread data for n threads max.

void TGeoVGXtru::CreateThreadData(Int_t nthreads)
{
   std::lock_guard<std::mutex> guard(fMutex);
   fThreadData.resize(nthreads);
   fThreadSize = nthreads;
   for (Int_t tid = 0; tid < nthreads; tid++) {
      if (fThreadData[tid] == nullptr) {
         fThreadData[tid] = new ThreadData_t;
         ThreadData_t &td = *fThreadData[tid];
         td.fXc = new Double_t[fNvert];
         td.fYc = new Double_t[fNvert];
         memcpy(td.fXc, fX, fNvert * sizeof(Double_t));
         memcpy(td.fYc, fY, fNvert * sizeof(Double_t));
         td.fPoly = new TGeoPolygon(fNvert);
         td.fPoly->SetXY(td.fXc, td.fYc); // initialize with current coordinates
         td.fPoly->FinishPolygon();
         if (tid == 0 && td.fPoly->IsIllegalCheck()) {
            Error("DefinePolygon", "Shape %s of type XTRU has an illegal polygon.", GetName());
         }
      }
   }
}

// ////////////////////////////////////////////////////////////////////////////////
// /// Set current z-plane.

// void TGeoVGXtru::SetIz(Int_t iz)
// {
//    GetThreadData().fIz = iz;
// }
// ////////////////////////////////////////////////////////////////////////////////
// /// Set current segment.

// void TGeoVGXtru::SetSeg(Int_t iseg)
// {
//    GetThreadData().fSeg = iseg;
// }

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoVGXtru::TGeoVGXtru()
   : Base_t(""),
     fNvert(0),
     fNz(0),
     fZcurrent(0.),
     fX(nullptr),
     fY(nullptr),
     fZ(nullptr),
     fScale(nullptr),
     fX0(nullptr),
     fY0(nullptr)
{
   SetShapeBit(TGeoShape::kGeoXtru);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGXtru::TGeoVGXtru(Int_t nz)
   : Base_t(""),
     fNvert(0),
     fNz(nz),
     fZcurrent(0.),
     fX(nullptr),
     fY(nullptr),
     fZ(new Double_t[nz]),
     fScale(new Double_t[nz]),
     fX0(new Double_t[nz]),
     fY0(new Double_t[nz])
{
   SetShapeBit(TGeoShape::kGeoXtru);
   if (nz < 2) {
      Error("ctor", "Cannot create TGeoVGXtru %s with less than 2 Z planes", GetName());
      SetShapeBit(TGeoShape::kGeoBad);
      return;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor in GEANT3 style
///  - param[0] = nz  // number of z planes
///
///  - param[1] = z1  // Z position of first plane
///  - param[2] = x1  // X position of first plane
///  - param[3] = y1  // Y position of first plane
///  - param[4] = scale1  // scale factor for first plane
/// ...
///  - param[4*(nz-1]+1] = zn
///  - param[4*(nz-1)+2] = xn
///  - param[4*(nz-1)+3] = yn
///  - param[4*(nz-1)+4] = scalen

TGeoVGXtru::TGeoVGXtru(Double_t *param)
   : Base_t(""),
     fNvert(0),
     fNz(0),
     fZcurrent(0.),
     fX(nullptr),
     fY(nullptr),
     fZ(nullptr),
     fScale(nullptr),
     fX0(nullptr),
     fY0(nullptr)
{
   SetShapeBit(TGeoShape::kGeoXtru);
   SetDimensions(param);
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGXtru::~TGeoVGXtru()
{
   if (fX) {
      delete[] fX;
      fX = nullptr;
   }
   if (fY) {
      delete[] fY;
      fY = nullptr;
   }
   if (fZ) {
      delete[] fZ;
      fZ = nullptr;
   }
   if (fScale) {
      delete[] fScale;
      fScale = nullptr;
   }
   if (fX0) {
      delete[] fX0;
      fX0 = nullptr;
   }
   if (fY0) {
      delete[] fY0;
      fY0 = nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the pcon

void TGeoVGXtru::ComputeBBox()
{
   ThreadData_t &td = GetThreadData();
   if (!fX || !fZ || !fNvert) {
      Error("ComputeBBox", "In shape %s polygon not defined", GetName());
      SetShapeBit(TGeoShape::kGeoBad);
      return;
   }
   Double_t zmin = fZ[0];
   Double_t zmax = fZ[fNz - 1];
   Double_t xmin = TGeoShape::Big();
   Double_t xmax = -TGeoShape::Big();
   Double_t ymin = TGeoShape::Big();
   Double_t ymax = -TGeoShape::Big();
   for (Int_t i = 0; i < fNz; i++) {
      SetCurrentVertices(fX0[i], fY0[i], fScale[i]);
      for (Int_t j = 0; j < fNvert; j++) {
         if (td.fXc[j] < xmin)
            xmin = td.fXc[j];
         if (td.fXc[j] > xmax)
            xmax = td.fXc[j];
         if (td.fYc[j] < ymin)
            ymin = td.fYc[j];
         if (td.fYc[j] > ymax)
            ymax = td.fYc[j];
      }
   }

   Double_t dx, dy, dz;
   Double_t origin[3];
   dx = 0.5 * (xmax - xmin);
   dy = 0.5 * (ymax - ymin);
   dz = 0.5 * (zmax - zmin);
   origin[0] = 0.5 * (xmin + xmax);
   origin[1] = 0.5 * (ymin + ymax);
   origin[2] = 0.5 * (zmin + zmax);
   fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);
}

////////////////////////////////////////////////////////////////////////////////
/// Divide this shape along one axis.

TGeoVolume *TGeoVGXtru::Divide(TGeoVolume *voldiv, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Division of a general trapezoid not implemented");
   return voldiv;
}

// ////////////////////////////////////////////////////////////////////////////////
// /// compute closest distance from point px,py to each corner

// Int_t TGeoVGXtru::DistancetoPrimitive(Int_t px, Int_t py)
// {
//    const Int_t numPoints = fNvert * fNz;
//    return ShapeDistancetoPrimitive(numPoints, px, py);
// }

////////////////////////////////////////////////////////////////////////////////
/// Draw the section polygon.

void TGeoVGXtru::DrawPolygon(Option_t *option)
{
   ThreadData_t &td = GetThreadData();
   if (td.fPoly)
      td.fPoly->Draw(option);
}

////////////////////////////////////////////////////////////////////////////////
/// Creates the polygon representing the blueprint of any Xtru section.
///  - nvert     = number of vertices >2
///  - xv[nvert] = array of X vertex positions
///  - yv[nvert] = array of Y vertex positions
///
/// *NOTE* should be called before DefineSection or ctor with 'param'

Bool_t TGeoVGXtru::DefinePolygon(Int_t nvert, const Double_t *xv, const Double_t *yv)
{
   if (nvert < 3) {
      Error("DefinePolygon", "In shape %s cannot create polygon with less than 3 vertices", GetName());
      SetShapeBit(TGeoShape::kGeoBad);
      return kFALSE;
   }
   for (Int_t i = 0; i < nvert - 1; i++) {
      for (Int_t j = i + 1; j < nvert; j++) {
         if (TMath::Abs(xv[i] - xv[j]) < TGeoShape::Tolerance() && TMath::Abs(yv[i] - yv[j]) < TGeoShape::Tolerance()) {
            Error("DefinePolygon", "In shape %s 2 vertices cannot be identical", GetName());
            SetShapeBit(TGeoShape::kGeoBad);
            //             return kFALSE;
         }
      }
   }
   fNvert = nvert;
   if (fX)
      delete[] fX;
   fX = new Double_t[nvert];
   if (fY)
      delete[] fY;
   fY = new Double_t[nvert];
   memcpy(fX, xv, nvert * sizeof(Double_t));
   memcpy(fY, yv, nvert * sizeof(Double_t));

   ClearThreadData();

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// defines z position of a section plane, rmin and rmax at this z.

void TGeoVGXtru::DefineSection(Int_t snum, Double_t z, Double_t x0, Double_t y0, Double_t scale)
{
   if ((snum < 0) || (snum >= fNz))
      return;
   fZ[snum] = z;
   fX0[snum] = x0;
   fY0[snum] = y0;
   fScale[snum] = scale;
   if (snum) {
      if (fZ[snum] < fZ[snum - 1]) {
         Warning("DefineSection",
                 "In shape: %s, Z position of section "
                 "%i, z=%e, not in increasing order, %i, z=%e",
                 GetName(), snum, fZ[snum], snum - 1, fZ[snum - 1]);
         return;
      }
   }
   if (snum == (fNz - 1)) {
      unsigned int nVertices = fNvert;
      unsigned int nSections = fNz;

      auto vertices = new vecgeom::XtruVertex2[nVertices];
      auto sections = new vecgeom::XtruSection[nSections];

      for (unsigned int i = 0; i < nVertices; ++i) {
         vertices[i].x = fX[i];
         vertices[i].y = fY[i];
      }
      for (unsigned int i = 0; i < nSections; ++i) {
         sections[i].fOrigin.Set(fX0[i],
                                 fY0[i],
                                 fZ[i]);
         sections[i].fScale = fScale[i];
      }
      Base_t::Initialize(nVertices, vertices, nSections, sections);
      ComputeBBox();
      if (TestShapeBit(TGeoShape::kGeoBad))
         InspectShape();
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Return the Z coordinate for segment ipl.

Double_t TGeoVGXtru::GetZ(Int_t ipl) const
{
   if (ipl < 0 || ipl > (fNz - 1)) {
      Error("GetZ", "In shape %s, ipl=%i out of range (0,%i)", GetName(), ipl, fNz - 1);
      return 0.;
   }
   return fZ[ipl];
}

////////////////////////////////////////////////////////////////////////////////
/// Print actual Xtru parameters.

void TGeoVGXtru::InspectShape() const
{
   printf("*** Shape %s: TGeoVGXtru ***\n", GetName());
   printf("    Nz    = %i\n", fNz);
   printf("    List of (x,y) of polygon vertices:\n");
   for (Int_t ivert = 0; ivert < fNvert; ivert++)
      printf("    x = %11.5f  y = %11.5f\n", fX[ivert], fY[ivert]);
   for (Int_t ipl = 0; ipl < fNz; ipl++)
      printf("     plane %i: z=%11.5f x0=%11.5f y0=%11.5f scale=%11.5f\n", ipl, fZ[ipl], fX0[ipl], fY0[ipl],
             fScale[ipl]);
   printf(" Bounding box:\n");
   fBoundingBox.InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGXtru::MakeBuffer3D() const
{
   Int_t nz = GetNz();
   Int_t nvert = GetNvert();
   Int_t nbPnts = nz * nvert;
   Int_t nbSegs = nvert * (2 * nz - 1);
   Int_t nbPols = nvert * (nz - 1) + 2;

   TBuffer3D *buff = new TBuffer3D(TBuffer3DTypes::kGeneric, nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols,
                                   6 * (nbPols - 2) + 2 * (2 + nvert));
   if (buff) {
      SetPoints(buff->fPnts);
      SetSegsAndPols(*buff);
   }

   return buff;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill TBuffer3D structure for segments and polygons.

void TGeoVGXtru::SetSegsAndPols(TBuffer3D &buff) const
{
   Int_t nz = GetNz();
   Int_t nvert = GetNvert();
   Int_t c = GetBasicColor();

   Int_t i, j;
   Int_t indx = 0, indx2, k;
   for (i = 0; i < nz; i++) {
      // loop Z planes
      indx2 = i * nvert;
      // loop polygon segments
      for (j = 0; j < nvert; j++) {
         k = (j + 1) % nvert;
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j;
         buff.fSegs[indx++] = indx2 + k;
      }
   } // total: nz*nvert polygon segments
   for (i = 0; i < nz - 1; i++) {
      // loop Z planes
      indx2 = i * nvert;
      // loop polygon segments
      for (j = 0; j < nvert; j++) {
         k = j + nvert;
         buff.fSegs[indx++] = c;
         buff.fSegs[indx++] = indx2 + j;
         buff.fSegs[indx++] = indx2 + k;
      }
   } // total (nz-1)*nvert lateral segments

   indx = 0;

   // fill lateral polygons
   for (i = 0; i < nz - 1; i++) {
      indx2 = i * nvert;
      for (j = 0; j < nvert; j++) {
         k = (j + 1) % nvert;
         buff.fPols[indx++] = c + j % 3;
         buff.fPols[indx++] = 4;
         buff.fPols[indx++] = indx2 + j;
         buff.fPols[indx++] = nz * nvert + indx2 + k;
         buff.fPols[indx++] = indx2 + nvert + j;
         buff.fPols[indx++] = nz * nvert + indx2 + j;
      }
   } // total (nz-1)*nvert polys
   buff.fPols[indx++] = c + 2;
   buff.fPols[indx++] = nvert;
   indx2 = 0;
   for (j = nvert - 1; j >= 0; --j) {
      buff.fPols[indx++] = indx2 + j;
   }

   buff.fPols[indx++] = c;
   buff.fPols[indx++] = nvert;
   indx2 = (nz - 1) * nvert;

   for (j = 0; j < nvert; j++) {
      buff.fPols[indx++] = indx2 + j;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGXtru::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   auto " << GetPointerName() << " = new TGeoVGXtru(" << fNz << ");" << std::endl;
   out << "   " << GetPointerName() << "->SetName(\"" << GetName() << "\");" << std::endl;
   for (Int_t i = 0; i < fNvert; i++) {
      out << "   xvert[" << i << "] = " << fX[i] << ";   yvert[" << i << "] = " << fY[i] << ";" << std::endl;
   }
   out << "   " << GetPointerName() << "->DefinePolygon(" << fNvert << ", xvert, yvert);" << std::endl;
   for (Int_t i = 0; i < fNz; i++)
      out << "   " << GetPointerName() << "->DefineSection(" << i << ", " << fZ[i] << ", " << fX0[i] << ", " << fY0[i]
          << ", " << fScale[i] << ");" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Recompute current section vertices for a given Z position within range of section iz.

void TGeoVGXtru::SetCurrentZ(Double_t z, Int_t iz)
{
   Double_t x0, y0, scale, a, b;
   Int_t ind1, ind2;
   ind1 = iz;
   ind2 = iz + 1;
   Double_t invdz = 1. / (fZ[ind2] - fZ[ind1]);
   a = (fX0[ind1] * fZ[ind2] - fX0[ind2] * fZ[ind1]) * invdz;
   b = (fX0[ind2] - fX0[ind1]) * invdz;
   x0 = a + b * z;
   a = (fY0[ind1] * fZ[ind2] - fY0[ind2] * fZ[ind1]) * invdz;
   b = (fY0[ind2] - fY0[ind1]) * invdz;
   y0 = a + b * z;
   a = (fScale[ind1] * fZ[ind2] - fScale[ind2] * fZ[ind1]) * invdz;
   b = (fScale[ind2] - fScale[ind1]) * invdz;
   scale = a + b * z;
   SetCurrentVertices(x0, y0, scale);
}

////////////////////////////////////////////////////////////////////////////////
/// Set current vertex coordinates according X0, Y0 and SCALE.

void TGeoVGXtru::SetCurrentVertices(Double_t x0, Double_t y0, Double_t scale)
{
   ThreadData_t &td = GetThreadData();
   for (Int_t i = 0; i < fNvert; i++) {
      td.fXc[i] = scale * fX[i] + x0;
      td.fYc[i] = scale * fY[i] + y0;
   }
}

////////////////////////////////////////////////////////////////////////////////
///  - param[0] = nz  // number of z planes
///
///  - param[1] = z1  // Z position of first plane
///  - param[2] = x1  // X position of first plane
///  - param[3] = y1  // Y position of first plane
///  - param[4] = scale1  // scale factor for first plane
/// ...
///  - param[4*(nz-1]+1] = zn
///  - param[4*(nz-1)+2] = xn
///  - param[4*(nz-1)+3] = yn
///  - param[4*(nz-1)+4] = scalen

void TGeoVGXtru::SetDimensions(Double_t *param)
{
   fNz = (Int_t)param[0];
   if (fNz < 2) {
      Error("SetDimensions", "Cannot create TGeoVGXtru %s with less than 2 Z planes", GetName());
      SetShapeBit(TGeoShape::kGeoBad);
      return;
   }
   if (fZ)
      delete[] fZ;
   if (fScale)
      delete[] fScale;
   if (fX0)
      delete[] fX0;
   if (fY0)
      delete[] fY0;
   fZ = new Double_t[fNz];
   fScale = new Double_t[fNz];
   fX0 = new Double_t[fNz];
   fY0 = new Double_t[fNz];

   for (Int_t i = 0; i < fNz; i++)
      DefineSection(i, param[1 + 4 * i], param[2 + 4 * i], param[3 + 4 * i], param[4 + 4 * i]);
}

////////////////////////////////////////////////////////////////////////////////
/// create polycone mesh points

void TGeoVGXtru::SetPoints(Double_t *points) const
{
   ThreadData_t &td = GetThreadData();
   Int_t i, j;
   Int_t indx = 0;
   TGeoVGXtru *xtru = (TGeoVGXtru *)this;
   if (points) {
      for (i = 0; i < fNz; i++) {
         xtru->SetCurrentVertices(fX0[i], fY0[i], fScale[i]);
         if (td.fPoly->IsClockwise()) {
            for (j = 0; j < fNvert; j++) {
               points[indx++] = td.fXc[j];
               points[indx++] = td.fYc[j];
               points[indx++] = fZ[i];
            }
         } else {
            for (j = 0; j < fNvert; j++) {
               points[indx++] = td.fXc[fNvert - 1 - j];
               points[indx++] = td.fYc[fNvert - 1 - j];
               points[indx++] = fZ[i];
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// create polycone mesh points

void TGeoVGXtru::SetPoints(Float_t *points) const
{
   ThreadData_t &td = GetThreadData();
   Int_t i, j;
   Int_t indx = 0;
   TGeoVGXtru *xtru = (TGeoVGXtru *)this;
   if (points) {
      for (i = 0; i < fNz; i++) {
         xtru->SetCurrentVertices(fX0[i], fY0[i], fScale[i]);
         if (td.fPoly->IsClockwise()) {
            for (j = 0; j < fNvert; j++) {
               points[indx++] = td.fXc[j];
               points[indx++] = td.fYc[j];
               points[indx++] = fZ[i];
            }
         } else {
            for (j = 0; j < fNvert; j++) {
               points[indx++] = td.fXc[fNvert - 1 - j];
               points[indx++] = td.fYc[fNvert - 1 - j];
               points[indx++] = fZ[i];
            }
         }
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGXtru::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   Int_t nz = GetNz();
   Int_t nv = GetNvert();
   nvert = nz * nv;
   nsegs = nv * (2 * nz - 1);
   npols = nv * (nz - 1) + 2;
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGXtru::GetNmeshVertices() const
{
   Int_t numPoints = fNz * fNvert;
   return numPoints;
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGXtru::Sizeof3D() const {}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGXtru::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
{
   static TBuffer3D buffer(TBuffer3DTypes::kGeneric);

   fBoundingBox.FillBuffer3D(buffer, reqSections, localFrame);

   if (reqSections & TBuffer3D::kRawSizes) {
      Int_t nz = GetNz();
      Int_t nvert = GetNvert();
      Int_t nbPnts = nz * nvert;
      Int_t nbSegs = nvert * (2 * nz - 1);
      Int_t nbPols = nvert * (nz - 1) + 2;
      if (buffer.SetRawSizes(nbPnts, 3 * nbPnts, nbSegs, 3 * nbSegs, nbPols, 6 * (nbPols - 2) + 2 * (2 + nvert))) {
         buffer.SetSectionsValid(TBuffer3D::kRawSizes);
      }
   }
   // TODO: Push down to TGeoShape?
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
