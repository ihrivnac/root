// @(#)root/geom:$Id$
// Author: Mihaela Gheata   05/06/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoEltu.h"
#include "TGeoVGEltu.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGEltu);

//////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGEltu::TGeoVGEltu()
 : Base_t("", 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoEltu);
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying X and Y semiaxis length

TGeoVGEltu::TGeoVGEltu(Double_t a, Double_t b, Double_t dz)
  : Base_t("", a, b, dz)
{
   SetShapeBit(TGeoShape::kGeoEltu);
   SetEltuDimensions(a, b, dz);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying X and Y semiaxis length

TGeoVGEltu::TGeoVGEltu(const char *name, Double_t a, Double_t b, Double_t dz)
 : Base_t(name, a, b, dz)
{
   SetName(name);
   SetShapeBit(TGeoShape::kGeoEltu);
   SetEltuDimensions(a, b, dz);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor specifying minimum and maximum radius
/// param[0] =  A
/// param[1] =  B
/// param[2] = dz

TGeoVGEltu::TGeoVGEltu(Double_t *param)
  : Base_t("", param[0], param[1], param[2])
{
   SetShapeBit(TGeoShape::kGeoEltu);
   SetDimensions(param);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGEltu::~TGeoVGEltu() {}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the tube

void TGeoVGEltu::ComputeBBox()
{
   fDX = Base_t::GetDx();
   fDY = Base_t::GetDy();
   fDZ = Base_t::GetDz();
}

////////////////////////////////////////////////////////////////////////////////
/// Divide the shape along one axis.

TGeoVolume *TGeoVGEltu::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/, Int_t /*ndiv*/,
                             Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Elliptical tubes divisions not implemented");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGEltu::GetBoundingCylinder(Double_t *param) const
{
   param[0] = 0.;                       // Rmin
   param[1] = TMath::Max(Base_t::GetDx(), Base_t::GetDy()); // Rmax
   param[1] *= param[1];
   param[2] = 0.;   // Phi1
   param[3] = 360.; // Phi2
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGEltu::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (!mother->TestShapeBit(kGeoEltu)) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t a, b, dz;
   a = Base_t::GetDx();
   b = Base_t::GetDy();
   dz = Base_t::GetDz();
   if (Base_t::GetDz() < 0)
      dz = ((TGeoVGEltu *)mother)->GetDz();
   if (Base_t::GetDx() < 0)
      a = ((TGeoVGEltu *)mother)->GetA();
   if (Base_t::GetDy() < 0)
      a = ((TGeoVGEltu *)mother)->GetB();

   return (new TGeoVGEltu(a, b, dz));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGEltu::InspectShape() const
{
   printf("*** Shape %s: TGeoVGEltu ***\n", GetName());
   printf("    A    = %11.5f\n", Base_t::GetDx());
   printf("    B    = %11.5f\n", Base_t::GetDy());
   printf("    dz   = %11.5f\n", Base_t::GetDz());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGEltu::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   a  = " << Base_t::GetDx() << ";" << std::endl;
   out << "   b  = " << Base_t::GetDy() << ";" << std::endl;
   out << "   dz = " << Base_t::GetDz() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGEltu(\"" << GetName() << "\",a,b,dz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions of the elliptical tube.

void TGeoVGEltu::SetEltuDimensions(Double_t a, Double_t b, Double_t dz)
{
   if ((a <= 0) || (b < 0) || (dz < 0)) {
      SetShapeBit(kGeoRunTimeShape);
   }
   Base_t::SetDx(a);
   Base_t::SetDy(b);
   Base_t::SetDz(dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Set shape dimensions starting from an array.

void TGeoVGEltu::SetDimensions(Double_t *param)
{
   Double_t a = param[0];
   Double_t b = param[1];
   Double_t dz = param[2];
   SetEltuDimensions(a, b, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Create elliptical tube mesh points

void TGeoVGEltu::SetPoints(Double_t *points) const
{
   Double_t dz;
   Int_t j, n;

   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   Double_t cph, sph;
   dz = Base_t::GetDz();

   Int_t indx = 0;
   Double_t r2, r;
   Double_t a2 = Base_t::GetDx() * Base_t::GetDx();
   Double_t b2 = Base_t::GetDy() * Base_t::GetDy();

   if (points) {
      for (j = 0; j < n; j++) {
         points[indx + 6 * n] = points[indx] = 0;
         indx++;
         points[indx + 6 * n] = points[indx] = 0;
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         sph = TMath::Sin(phi);
         cph = TMath::Cos(phi);
         r2 = (a2 * b2) / (b2 + (a2 - b2) * sph * sph);
         r = TMath::Sqrt(r2);
         points[indx + 6 * n] = points[indx] = r * cph;
         indx++;
         points[indx + 6 * n] = points[indx] = r * sph;
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGEltu::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   // TGeoTube::GetMeshNumbers(nvert, nsegs, npols);
   // copied from TGeoTube
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
/// Returns the number of vertices on the mesh.

Int_t TGeoVGEltu::GetNmeshVertices() const
{
   // return TGeoTube::GetNmeshVertices();
   // copied from TGeoTube
   Int_t n = gGeoManager->GetNsegments() + 1;
   Int_t numPoints = n * 4;
   return numPoints;

}

////////////////////////////////////////////////////////////////////////////////
/// Create elliptical tube mesh points

void TGeoVGEltu::SetPoints(Float_t *points) const
{
   Double_t dz;
   Int_t j, n;

   n = gGeoManager->GetNsegments();
   Double_t dphi = 360. / n;
   Double_t phi = 0;
   Double_t cph, sph;
   dz = Base_t::GetDz();

   Int_t indx = 0;
   Double_t r2, r;
   Double_t a2 = Base_t::GetDx() * Base_t::GetDx();
   Double_t b2 = Base_t::GetDy() * Base_t::GetDy();

   if (points) {
      for (j = 0; j < n; j++) {
         points[indx + 6 * n] = points[indx] = 0;
         indx++;
         points[indx + 6 * n] = points[indx] = 0;
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
      for (j = 0; j < n; j++) {
         phi = j * dphi * TMath::DegToRad();
         sph = TMath::Sin(phi);
         cph = TMath::Cos(phi);
         r2 = (a2 * b2) / (b2 + (a2 - b2) * sph * sph);
         r = TMath::Sqrt(r2);
         points[indx + 6 * n] = points[indx] = r * cph;
         indx++;
         points[indx + 6 * n] = points[indx] = r * sph;
         indx++;
         points[indx + 6 * n] = dz;
         points[indx] = -dz;
         indx++;
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills a static 3D buffer and returns a reference.

const TBuffer3D &TGeoVGEltu::GetBuffer3D(Int_t reqSections, Bool_t localFrame) const
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
