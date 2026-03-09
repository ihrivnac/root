// @(#)root/geom:$Id$// Author: Andrei Gheata   24/10/01

// Contains() and DistFromOutside/Out() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoBBox.h"
#include "TGeoVGBBox.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include <iostream>

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TVirtualGeoPainter.h"
#include "TGeoBaseBox.h"
#include "TBuffer3D.h"
#include "TBuffer3DTypes.h"
#include "TMath.h"
#include "TRandom.h"

// ClassImp(TGeoVGBBox);

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGBBox::TGeoVGBBox()
  : Base_t("", 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoBox);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor where half-lengths are provided.

TGeoVGBBox::TGeoVGBBox(Double_t dx, Double_t dy, Double_t dz)
  : Base_t("", dx, dy, dz)
{
   SetShapeBit(TGeoShape::kGeoBox);
   SetBoxDimensions(dx, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with shape name.

TGeoVGBBox::TGeoVGBBox(const char *name, Double_t dx, Double_t dy, Double_t dz)
  : Base_t(name, dx, dy, dz)
{
   SetShapeBit(TGeoShape::kGeoBox);
   SetBoxDimensions(dx, dy, dz);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor based on the array of parameters.
///  - param[0] - half-length in x
///  - param[1] - half-length in y
///  - param[2] - half-length in z

TGeoVGBBox::TGeoVGBBox(Double_t *param)
  : Base_t("", 0., 0., 0.)
{
   SetShapeBit(TGeoShape::kGeoBox);
   SetDimensions(param);
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

TGeoVGBBox::~TGeoVGBBox() {}


////////////////////////////////////////////////////////////////////////////////
/// Divide this box shape belonging to volume "voldiv" into ndiv equal volumes
/// called divname, from start position with the given step. Returns pointer
/// to created division cell volume. In case a wrong division axis is supplied,
/// returns pointer to volume to be divided.

TGeoVolume *
TGeoVGBBox::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Double_t end = start + ndiv * step;
   switch (iaxis) {
   case 1: //--- divide on X
      shape = new TGeoVGBBox(step / 2., Base_t::y(), Base_t::z());
      finder = new TGeoPatternX(voldiv, ndiv, start, end);
      opt = "X";
      break;
   case 2: //--- divide on Y
      shape = new TGeoVGBBox(Base_t::x(), step / 2., Base_t::z());
      finder = new TGeoPatternY(voldiv, ndiv, start, end);
      opt = "Y";
      break;
   case 3: //--- divide on Z
      shape = new TGeoVGBBox(Base_t::x(), Base_t::y(), step / 2.);
      finder = new TGeoPatternZ(voldiv, ndiv, start, end);
      opt = "Z";
      break;
   default: Error("Divide", "Wrong axis type for division"); return nullptr;
   }
   vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
   vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
   vmulti->AddVolume(vol);
   voldiv->SetFinder(finder);
   finder->SetDivIndex(voldiv->GetNdaughters());
   for (Int_t ic = 0; ic < ndiv; ic++) {
      voldiv->AddNodeOffset(vol, ic, start + step / 2. + ic * step, opt.Data());
      ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
   }
   return vmulti;
}

////////////////////////////////////////////////////////////////////////////////
/// Static method to check if point[3] is located inside a box of having dx, dy, dz
/// as half-lengths.

Bool_t TGeoVGBBox::Contains(const Double_t *point, Double_t dx, Double_t dy, Double_t dz, const Double_t *origin)
{
   if (TMath::Abs(point[2] - origin[2]) > dz)
      return kFALSE;
   if (TMath::Abs(point[0] - origin[0]) > dx)
      return kFALSE;
   if (TMath::Abs(point[1] - origin[1]) > dy)
      return kFALSE;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// Static: Compute distance from outside point to surface of the box.
/// Boundary safe algorithm.

Double_t TGeoVGBBox::DistFromOutside(const Double_t *point, const Double_t *dir, Double_t dx, Double_t dy, Double_t dz,
                                   const Double_t *origin, Double_t stepmax)
{
   Bool_t in = kTRUE;
   Double_t saf[3];
   Double_t par[3];
   Double_t newpt[3];
   Int_t i, j;
   for (i = 0; i < 3; i++)
      newpt[i] = point[i] - origin[i];
   par[0] = dx;
   par[1] = dy;
   par[2] = dz;
   for (i = 0; i < 3; i++) {
      saf[i] = TMath::Abs(newpt[i]) - par[i];
      if (saf[i] >= stepmax)
         return TGeoShape::Big();
      if (in && saf[i] > 0)
         in = kFALSE;
   }
   // In case point is inside return ZERO
   if (in)
      return 0.0;
   Double_t coord, snxt = TGeoShape::Big();
   Int_t ibreak = 0;
   for (i = 0; i < 3; i++) {
      if (saf[i] < 0)
         continue;
      if (newpt[i] * dir[i] >= 0)
         continue;
      snxt = saf[i] / TMath::Abs(dir[i]);
      ibreak = 0;
      for (j = 0; j < 3; j++) {
         if (j == i)
            continue;
         coord = newpt[j] + snxt * dir[j];
         if (TMath::Abs(coord) > par[j]) {
            ibreak = 1;
            break;
         }
      }
      if (!ibreak)
         return snxt;
   }
   return TGeoShape::Big();
}

////////////////////////////////////////////////////////////////////////////////
/// In case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGBBox::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   Double_t dx, dy, dz;
   TGeoBaseBox baseBox(Base_t::x(), Base_t::y(), Base_t::z());
   Int_t ierr = mother->GetFittingBox(&baseBox, mat, dx, dy, dz);
   if (ierr) {
      Error("GetMakeRuntimeShape", "cannot fit this to mother");
      return nullptr;
   }
   return (new TGeoVGBBox(dx, dy, dz));
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGBBox::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   nvert = 8;
   nsegs = 12;
   npols = 6;
}

////////////////////////////////////////////////////////////////////////////////
/// Prints shape parameters

void TGeoVGBBox::InspectShape() const
{
   printf("*** Shape %s: TGeoVGBBox ***\n", GetName());
   printf("    dX = %11.5f\n", Base_t::x());
   printf("    dY = %11.5f\n", Base_t::y());
   printf("    dZ = %11.5f\n", Base_t::z());
}

////////////////////////////////////////////////////////////////////////////////
/// Creates a TBuffer3D describing *this* shape.
/// Coordinates are in local reference frame.

TBuffer3D *TGeoVGBBox::MakeBuffer3D() const
{
   TBuffer3D *buff = new TBuffer3D(TBuffer3DTypes::kGeneric, 8, 24, 12, 36, 6, 36);
   if (buff) {
      SetPoints(buff->fPnts);
      SetSegsAndPols(*buff);
   }

   return buff;
}

///
////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGBBox::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dx = " << Base_t::x() << ";" << std::endl;
   out << "   dy = " << Base_t::y() << ";" << std::endl;
   out << "   dz = " << Base_t::z() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGBBox(\"" << GetName() << "\", dx,dy,dz);" << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// Set parameters of the box.

void TGeoVGBBox::SetBoxDimensions(Double_t dx, Double_t dy, Double_t dz)
{
   // Set data to VecGeom base
   SetX(dx);
   SetY(dy);
   SetZ(dz);

   // Set bounding box
   Double_t origin[3];
   memset(origin, 0, 3 * sizeof(Double_t));
   fBoundingBox.SetBoxDimensions(dx, dy, dz, origin);

   if (TMath::Abs(Base_t::x()) < TGeoShape::Tolerance() && TMath::Abs(Base_t::y()) < TGeoShape::Tolerance() &&
       TMath::Abs(Base_t::z()) < TGeoShape::Tolerance())
      return;
   if ((dx < 0) || (dy < 0) || (dz < 0))
      SetShapeBit(kGeoRunTimeShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Set dimensions based on the array of parameters
/// param[0] - half-length in x
/// param[1] - half-length in y
/// param[2] - half-length in z

void TGeoVGBBox::SetDimensions(Double_t *param)
{
   if (!param) {
      Error("SetDimensions", "null parameters");
      return;
   }
   SetBoxDimensions(param[0], param[1], param[2]);
   Double_t origin[3];
   memset(origin, 0, 3 * sizeof(Double_t));   
   fBoundingBox.SetBoxDimensions(param[0], param[1], param[2], origin);

   if (TMath::Abs(Base_t::x()) < TGeoShape::Tolerance() && TMath::Abs(Base_t::y()) < TGeoShape::Tolerance() &&
       TMath::Abs(Base_t::z()) < TGeoShape::Tolerance())
      return;
   if ((Base_t::x() < 0) || (Base_t::y() < 0) || (Base_t::z() < 0))
      SetShapeBit(kGeoRunTimeShape);
}

////////////////////////////////////////////////////////////////////////////////
/// Fill box vertices to an array.

void TGeoVGBBox::SetBoxPoints(Double_t *points) const
{
   Base_t::SetPoints(points);
}

#endif
