// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02
// TGeoVGTrd::Contains() and DistFromInside() implemented by Mihaela Gheata

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoVGTrd.h"

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVolume.h"
#include "TMath.h"

#include <iostream>

// ClassImp(TGeoVGTrd);

////////////////////////////////////////////////////////////////////////////////
/// dummy ctor

TGeoVGTrd::TGeoVGTrd()
  : Base_t("", 0., 0., 0., 0., 0.)
{
   SetShapeBit(kGeoTrd2);
}

////////////////////////////////////////////////////////////////////////////////
/// constructor (Trd1).

TGeoVGTrd::TGeoVGTrd(Double_t dx1, Double_t dx2, Double_t dy, Double_t dz)
  : Base_t("", dx1, dx2, dy, dz)
{
   SetShapeBit(kGeoTrd1);
   if ((dx1 < 0) || (dx2 < 0) || (dy < 0) || (dz < 0)) {
      SetShapeBit(kGeoRunTimeShape);
      printf("trd1 : dx1=%f, dx2=%f, dy=%f, dz=%f\n", dx1, dx2, dy, dz);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor  (Trd1).

TGeoVGTrd::TGeoVGTrd(const char *name, Double_t dx1, Double_t dx2, Double_t dy, Double_t dz)
  : Base_t(name, dx1, dx2, dy, dz)
{
   SetShapeBit(kGeoTrd1);
   if ((dx1 < 0) || (dx2 < 0) || (dy < 0) || (dz < 0)) {
      SetShapeBit(kGeoRunTimeShape);
      printf("trd1 : dx1=%f, dx2=%f, dy=%f, dz=%f\n", dx1, dx2, dy, dz);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor  (Trd2).

TGeoVGTrd::TGeoVGTrd(Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2, Double_t dz)
  : Base_t("", dx1, dx2, dy1, dy2, dz)
{
   SetShapeBit(kGeoTrd2);
   if ((dx1 < 0) || (dx2 < 0) || (dy1 < 0) || (dy2 < 0) || (dz < 0)) {
      SetShapeBit(kGeoRunTimeShape);
      printf("trd2 : dx1=%f, dx2=%f, dy1=%f, dy2=%f, dz=%f\n", dx1, dx2, dy1, dy2, dz);
   } else
      ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// constructor  (Trd2).

TGeoVGTrd::TGeoVGTrd(const char *name, Double_t dx1, Double_t dx2, Double_t dy1, Double_t dy2, Double_t dz)
  : Base_t(name, dx1, dx2, dy1, dy2, dz)
{
   SetShapeBit(kGeoTrd2);
   if ((dx1 < 0) || (dx2 < 0) || (dy1 < 0) || (dy2 < 0) || (dz < 0)) {
      SetShapeBit(kGeoRunTimeShape);
      printf("trd2 : dx1=%f, dx2=%f, dy1=%f, dy2=%f, dz=%f\n", dx1, dx2, dy1, dy2, dz);
   } else
      ComputeBBox();
}

// ////////////////////////////////////////////////////////////////////////////////
// /// ctor with an array of parameters
// ///  - param[0] = dx1
// ///  - param[1] = dx2
// ///  - param[2] = dy1
// ///  - param[3] = dy2
// ///  - param[4] = dz

// TGeoVGTrd::TGeoVGTrd(Double_t *param)
//   :  Base_t("", param[0], param[1], param[2], param[3], param[4])
// {
//    SetShapeBit(kGeoTrd2);
//    SetDimensions(param);
//    if ((dx1() < 0) || (dx2() < 0) || (dy1() < 0) || (dy2() < 0) || (dz() < 0))
//       SetShapeBit(kGeoRunTimeShape);
//    else
//       ComputeBBox();
// }

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGTrd::~TGeoVGTrd() {}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box for a trd2

void TGeoVGTrd::ComputeBBox()
{
   fDX = TMath::Max(dx1(), dx2());
   fDY = TMath::Max(dy1(), dy2());
   fDZ = dz();
   memset(fOrigin, 0, 3 * sizeof(Double_t));
}

////////////////////////////////////////////////////////////////////////////////
/// Get range of shape for a given axis.

Double_t TGeoVGTrd::GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const
{
   xlo = 0;
   xhi = 0;
   Double_t dx = 0;
   switch (iaxis) {
   case 3:
      xlo = -dz();
      xhi = dz();
      dx = xhi - xlo;
      return dx;
   }
   return dx;
}

// ////////////////////////////////////////////////////////////////////////////////
// /// get the most visible corner from outside point and the normals

// void TGeoVGTrd::GetVisibleCorner(const Double_t *point, Double_t *vertex, Double_t *normals) const
// {
//    Double_t fx = 0.5 * (fDx1 - fDx2) / fDz;
//    Double_t fy = 0.5 * (fDy1 - fDy2) / fDz;
//    Double_t calf = 1. / TMath::Sqrt(1.0 + fx * fx);
//    Double_t salf = calf * fx;
//    Double_t cbet = 1. / TMath::Sqrt(1.0 + fy * fy);
//    Double_t sbet = cbet * fy;
//    // check visibility of X,Y faces
//    Double_t distx = fDx1 - fx * (fDz + point[2]);
//    Double_t disty = fDy1 - fy * (fDz + point[2]);
//    memset(normals, 0, 9 * sizeof(Double_t));
//    TGeoVGTrd *trd2 = (TGeoVGTrd *)this;
//    if (point[0] > distx) {
//       // hi x face visible
//       trd2->SetShapeBit(kGeoVisX);
//       normals[0] = calf;
//       normals[2] = salf;
//    } else {
//       trd2->SetShapeBit(kGeoVisX, kFALSE);
//       normals[0] = -calf;
//       normals[2] = salf;
//    }
//    if (point[1] > disty) {
//       // hi y face visible
//       trd2->SetShapeBit(kGeoVisY);
//       normals[4] = cbet;
//       normals[5] = sbet;
//    } else {
//       trd2->SetShapeBit(kGeoVisY, kFALSE);
//       normals[4] = -cbet;
//       normals[5] = sbet;
//    }
//    if (point[2] > fDz) {
//       // hi z face visible
//       trd2->SetShapeBit(kGeoVisZ);
//       normals[8] = 1;
//    } else {
//       trd2->SetShapeBit(kGeoVisZ, kFALSE);
//       normals[8] = -1;
//    }
//    SetVertex(vertex);
// }

// ////////////////////////////////////////////////////////////////////////////////
// /// get the opposite corner of the intersected face

// void TGeoVGTrd::GetOppositeCorner(const Double_t * /*point*/, Int_t inorm, Double_t *vertex, Double_t *normals) const
// {
//    TGeoVGTrd *trd2 = (TGeoVGTrd *)this;
//    if (inorm != 0) {
//       // change x face
//       trd2->SetShapeBit(kGeoVisX, !TestShapeBit(kGeoVisX));
//       normals[0] = -normals[0];
//    }
//    if (inorm != 1) {
//       // change y face
//       trd2->SetShapeBit(kGeoVisY, !TestShapeBit(kGeoVisY));
//       normals[4] = -normals[4];
//    }
//    if (inorm != 2) {
//       // hi z face visible
//       trd2->SetShapeBit(kGeoVisZ, !TestShapeBit(kGeoVisZ));
//       normals[8] = -normals[8];
//    }
//    SetVertex(vertex);
// }

////////////////////////////////////////////////////////////////////////////////
/// Divide this trd2 shape belonging to volume "voldiv" into ndiv volumes
/// called divname, from start position with the given step. Only Z divisions
/// are supported. For Z divisions just return the pointer to the volume to be
/// divided. In case a wrong division axis is supplied, returns pointer to
/// volume that was divided.

TGeoVolume *
TGeoVGTrd::Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start, Double_t step)
{
   TGeoShape *shape;          //--- shape to be created
   TGeoVolume *vol;           //--- division volume to be created
   TGeoVolumeMulti *vmulti;   //--- generic divided volume
   TGeoPatternFinder *finder; //--- finder to be attached
   TString opt = "";          //--- option to be attached
   Double_t zmin, zmax, dx1n, dx2n, dy1n, dy2n;
   Int_t id;
   Double_t end = start + ndiv * step;
   switch (iaxis) {
   case 1: Warning("Divide", "dividing a Trd2 on X not implemented"); return nullptr;
   case 2: Warning("Divide", "dividing a Trd2 on Y not implemented"); return nullptr;
   case 3:
      finder = new TGeoPatternZ(voldiv, ndiv, start, end);
      vmulti = gGeoManager->MakeVolumeMulti(divname, voldiv->GetMedium());
      voldiv->SetFinder(finder);
      finder->SetDivIndex(voldiv->GetNdaughters());
      for (id = 0; id < ndiv; id++) {
         zmin = start + id * step;
         zmax = start + (id + 1) * step;
         dx1n = 0.5 * (dx1() * (dz() - zmin) + dx2() * (dz() + zmin)) / dz();
         dx2n = 0.5 * (dx1() * (dz() - zmax) + dx2() * (dz() + zmax)) / dz();
         dy1n = 0.5 * (dy1() * (dz() - zmin) + dy2() * (dz() + zmin)) / dz();
         dy2n = 0.5 * (dy1() * (dz() - zmax) + dy2() * (dz() + zmax)) / dz();
         shape = new TGeoVGTrd(dx1n, dx2n, dy1n, dy2n, step / 2.);
         vol = new TGeoVolume(divname, shape, voldiv->GetMedium());
         vmulti->AddVolume(vol);
         opt = "Z";
         voldiv->AddNodeOffset(vol, id, start + step / 2 + id * step, opt.Data());
         ((TGeoNodeOffset *)voldiv->GetNodes()->At(voldiv->GetNdaughters() - 1))->SetFinder(finder);
      }
      return vmulti;
   default: Error("Divide", "Wrong axis type for division"); return nullptr;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Fill vector param[4] with the bounding cylinder parameters. The order
/// is the following : Rmin, Rmax, Phi1, Phi2

void TGeoVGTrd::GetBoundingCylinder(Double_t *param) const
{
   TGeoBBox::GetBoundingCylinder(param);
}

////////////////////////////////////////////////////////////////////////////////
/// Return dy() if Trd is of TRd1 type
Double_t TGeoVGTrd::GetDy() const
{
  if (dy1() == dy2()) {
     return dy1();
  }
  else {
     Error("GetDy", "The shape is not of Trd1 type");
     return 0.;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// Fills real parameters of a positioned box inside this. Returns 0 if successful.

Int_t TGeoVGTrd::GetFittingBox(const TGeoBBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const
{
   dx = dy = dz = 0;
   if (mat->IsRotation()) {
      Error("GetFittingBox", "cannot handle parametrized rotated volumes");
      return 1; // ### rotation not accepted ###
   }
   //--> translate the origin of the parametrized box to the frame of this box.
   Double_t origin[3];
   mat->LocalToMaster(parambox->GetOrigin(), origin);
   if (!Contains(origin)) {
      Error("GetFittingBox", "wrong matrix - parametrized box is outside this");
      return 1; // ### wrong matrix ###
   }
   //--> now we have to get the valid range for all parametrized axis
   Double_t dd[3];
   dd[0] = parambox->GetDX();
   dd[1] = parambox->GetDY();
   dd[2] = parambox->GetDZ();
   //-> check if Z range is fixed
   if (dd[2] < 0) {
      dd[2] = TMath::Min(origin[2] + Base_t::dz(), Base_t::dz() - origin[2]);
      if (dd[2] < 0) {
         Error("GetFittingBox", "wrong matrix");
         return 1;
      }
   }
   if (dd[0] >= 0 && dd[1] >= 0) {
      dx = dd[0];
      dy = dd[1];
      dz = dd[2];
      return 0;
   }
   //-> check now range at Z = origin[2] +/- dd[2]
   Double_t fx = 0.5 * (dx1() - dx2()) / Base_t::dz();
   Double_t fy = 0.5 * (dy1() - dy2()) / Base_t::dz();
   Double_t dx0 = 0.5 * (dx1() + dx2());
   Double_t dy0 = 0.5 * (dy1() + dy2());
   Double_t z = origin[2] - dd[2];
   dd[0] = dx0 - fx * z - origin[0];
   dd[1] = dy0 - fy * z - origin[1];
   z = origin[2] + dd[2];
   dd[0] = TMath::Min(dd[0], dx0 - fx * z - origin[0]);
   dd[1] = TMath::Min(dd[1], dy0 - fy * z - origin[1]);
   if (dd[0] < 0 || dd[1] < 0) {
      Error("GetFittingBox", "wrong matrix");
      return 1;
   }
   dx = dd[0];
   dy = dd[1];
   dz = dd[2];
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// in case shape has some negative parameters, these has to be computed
/// in order to fit the mother

TGeoShape *TGeoVGTrd::GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix * /*mat*/) const
{
   if (!TestShapeBit(kGeoRunTimeShape))
      return nullptr;
   if (!mother->TestShapeBit(kGeoTrd2)) {
      Error("GetMakeRuntimeShape", "invalid mother");
      return nullptr;
   }
   Double_t dx1, dx2, dy1, dy2, dz;
   if (Base_t::dx1() < 0)
      dx1 = ((TGeoVGTrd *)mother)->GetDx1();
   else
      dx1 = Base_t::dx1();
   if (Base_t::dx2() < 0)
      dx2 = ((TGeoVGTrd *)mother)->GetDx2();
   else
      dx2 = Base_t::dx2();
   if (Base_t::dy1() < 0)
      dy1 = ((TGeoVGTrd *)mother)->GetDy1();
   else
      dy1 = Base_t::dy1();
   if (Base_t::dy2() < 0)
      dy2 = ((TGeoVGTrd *)mother)->GetDy2();
   else
      dy2 = Base_t::dy2();
   if (Base_t::dz() < 0)
      dz = ((TGeoVGTrd *)mother)->GetDz();
   else
      dz = Base_t::dz();

   return (new TGeoVGTrd(dx1, dx2, dy1, dy2, dz));
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGTrd::InspectShape() const
{
   printf("*** Shape %s: TGeoVGTrd ***\n", GetName());
   printf("    dx1 = %11.5f\n", dx1());
   printf("    dx2 = %11.5f\n", dx2());
   printf("    dy1 = %11.5f\n", dy1());
   printf("    dy2 = %11.5f\n", dy2());
   printf("    dz  = %11.5f\n", dz());
   printf(" Bounding box:\n");
   TGeoBBox::InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGTrd::SavePrimitive(std::ostream &out, Option_t * /*option*/ /*= ""*/)
{
   if (TObject::TestBit(kGeoSavePrimitive))
      return;
   out << "   // Shape: " << GetName() << " type: " << ClassName() << std::endl;
   out << "   dx1 = " << dx1() << ";" << std::endl;
   out << "   dx2 = " << dx2() << ";" << std::endl;
   out << "   dy1 = " << dy1() << ";" << std::endl;
   out << "   dy2 = " << dy2() << ";" << std::endl;
   out << "   dz  = " << dz() << ";" << std::endl;
   out << "   TGeoShape *" << GetPointerName() << " = new TGeoVGTrd(\"" << GetName() << "\", dx1,dx2,dy1,dy2,dz);"
       << std::endl;
   TObject::SetBit(TGeoShape::kGeoSavePrimitive);
}

////////////////////////////////////////////////////////////////////////////////
/// set arb8 params in one step :

void TGeoVGTrd::SetDimensions(Double_t *param)
{
   SetXHalfLength1(param[0]);
   SetXHalfLength2(param[1]);
   SetYHalfLength1(param[2]);
   SetYHalfLength2(param[3]);
   SetZHalfLength(param[4]);
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// create trd2 mesh points

void TGeoVGTrd::SetPoints(Double_t *points) const
{
   if (!points)
      return;
   points[0] = -dx1();
   points[1] = -dy1();
   points[2] = -dz();
   points[3] = -dx1();
   points[4] = dy1();
   points[5] = -dz();
   points[6] = dx1();
   points[7] = dy1();
   points[8] = -dz();
   points[9] = dx1();
   points[10] = -dy1();
   points[11] = -dz();
   points[12] = -dx2();
   points[13] = -dy2();
   points[14] = dz();
   points[15] = -dx2();
   points[16] = dy2();
   points[17] = dz();
   points[18] = dx2();
   points[19] = dy2();
   points[20] = dz();
   points[21] = dx2();
   points[22] = -dy2();
   points[23] = dz();
}

////////////////////////////////////////////////////////////////////////////////
/// create trd2 mesh points

void TGeoVGTrd::SetPoints(Float_t *points) const
{
   if (!points)
      return;
   points[0] = -dx1();
   points[1] = -dy1();
   points[2] = -dz();
   points[3] = -dx1();
   points[4] = dy1();
   points[5] = -dz();
   points[6] = dx1();
   points[7] = dy1();
   points[8] = -dz();
   points[9] = dx1();
   points[10] = -dy1();
   points[11] = -dz();
   points[12] = -dx2();
   points[13] = -dy2();
   points[14] = dz();
   points[15] = -dx2();
   points[16] = dy2();
   points[17] = dz();
   points[18] = dx2();
   points[19] = dy2();
   points[20] = dz();
   points[21] = dx2();
   points[22] = -dy2();
   points[23] = dz();
}

////////////////////////////////////////////////////////////////////////////////
/// set vertex of a corner according to visibility flags

void TGeoVGTrd::SetVertex(Double_t *vertex) const
{
   if (TestShapeBit(kGeoVisX)) {
      if (TestShapeBit(kGeoVisZ)) {
         vertex[0] = dx2();
         vertex[2] = dz();
         vertex[1] = (TestShapeBit(kGeoVisY)) ? dy2() : -dy2();
      } else {
         vertex[0] = dx1();
         vertex[2] = -dz();
         vertex[1] = (TestShapeBit(kGeoVisY)) ? dy1() : -dy1();
      }
   } else {
      if (TestShapeBit(kGeoVisZ)) {
         vertex[0] = -dx2();
         vertex[2] = dz();
         vertex[1] = (TestShapeBit(kGeoVisY)) ? dy2() : -dy2();
      } else {
         vertex[0] = -dx1();
         vertex[2] = -dz();
         vertex[1] = (TestShapeBit(kGeoVisY)) ? dy1() : -dy1();
      }
   }
}

////////////////////////////////////////////////////////////////////////////////
/// fill size of this 3-D object

void TGeoVGTrd::Sizeof3D() const
{
   TGeoBBox::Sizeof3D();
}

#endif
