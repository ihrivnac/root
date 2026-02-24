// @(#)root/geom:$Id$
// Author: Andrei Gheata   31/01/02

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#include "TGeoVGMultiUnion.h"

#include <iostream>
#include <vector>

#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoBoolNode.h"


// ClassImp(TGeoVGMultiUnion);

namespace {

/**
 * Converts a TGeoHMatrix to a TGeoCombiTrans.
 * @param matrix The input HMatrix containing rotation and translation.
 * @return A new TGeoCombiTrans object representing the same transformation.
 */
// TGeoCombiTrans* ConvertHMatrixToCombiTrans(const TGeoHMatrix& matrix) {
//     // 1. Extract the translation components
//     const Double_t* translation = matrix.GetTranslation();
    
//     // 2. Extract the rotation matrix components (3x3)
//     const Double_t* rotation = matrix.GetRotationMatrix();

//     // 3. Create the CombiTrans object
//     // We initialize it with the translation vector
//     TGeoCombiTrans* combi = new TGeoCombiTrans(translation[0], translation[1], translation[2], nullptr);

//     // 4. Create a TGeoRotation object and set it
//     TGeoRotation rot;
//     rot.SetMatrix(rotation);
//     combi->SetRotation(rot);

//     return combi;
// }

////////////////////////////////////////////////////////////////////////////////
/// Convert a TGeoMatrix to a TRansformation3D
/// From TGeoVGShape

vecgeom::cxx::Transformation3D* Convert(TGeoMatrix const *const geomatrix)
{
   Double_t const *const t = geomatrix->GetTranslation();
   Double_t const *const r = geomatrix->GetRotationMatrix();
   vecgeom::cxx::Transformation3D *const transformation =
      new vecgeom::cxx::Transformation3D(t[0], t[1], t[2], r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]);
   return transformation;
}

// typedef std::pair<TGeoShape*, TGeoHMatrix*> ShapeEntry;

/**
 * Creates a TGeoCompositeShape equivalent to a MultiUnion by nesting TGeoUnion nodes.
 * @param name The name of the resulting composite shape.
 * @param components A vector of shapes and their relative transformations.
 * @return A pointer to the newly created TGeoCompositeShape.
 */
TGeoCompositeShape* MakeComposite(const char *name, const std::vector<ShapeEntry>& components) 
{
    if (components.empty()) return nullptr;
    if (components.size() == 1) {
        // A composite of one shape is just the shape with its matrix
        return new TGeoCompositeShape(name, new TGeoUnion(components[0].first, nullptr, 
                                      components[0].second, nullptr));
    }

    // Initialize the left side of the tree with the first component
    TGeoShape* leftShape = components[0].first;
    TGeoMatrix* leftMat  = components[0].second;

    // Iteratively build the tree: Result = ((A + B) + C) ...
    for (size_t i = 1; i < components.size(); ++i) {
        TGeoShape* rightShape = components[i].first;
        TGeoMatrix* rightMat  = components[i].second;

        // Create a new Union node joining the current tree (left) with the next shape (right)
        TGeoUnion* boolNode = new TGeoUnion(leftShape, rightShape, leftMat, rightMat);
        
        // Wrap it in a temporary composite to use as the 'left' side for the next iteration
        // We use a temporary name for internal nodes
        leftShape = new TGeoCompositeShape(Form("%s_branch_%lu", name, i), boolNode);
        
        // The transformation is now baked into the BoolNode, so the next leftMat is Identity
        leftMat = gGeoIdentity; 
    }

    // The final 'leftShape' is the complete tree
    return (TGeoCompositeShape*)leftShape;
}

}

////////////////////////////////////////////////////////////////////////////////
/// Default constructor

TGeoVGMultiUnion::TGeoVGMultiUnion()
 : Base_t("")
{
   SetShapeBit(TGeoShape::kGeoComb);
}

////////////////////////////////////////////////////////////////////////////////
/// Constructor with a components list

TGeoVGMultiUnion::TGeoVGMultiUnion(const char *name, const std::vector<ShapeEntry>& components)
  : Base_t(name)
{
   // Vecgeom MultiUnion
   for (auto& item : components) {
      // cast TGeoShape to VUnplacedVolume
      auto unplacedVolume = dynamic_cast<const VUnplacedVolume*>(item.first);
      if (unplacedVolume == nullptr) {
         std::cerr << "Component shape is not VecGeom type !!!" << std::endl;
      }
      AddNode(unplacedVolume, *Convert(item.second));
      delete item.second;
   }
   Close();

   // Composite shape
   fCompositeShape = MakeComposite(name, components);

   // original composite shape for non navigation functions
   ComputeBBox();
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

TGeoVGMultiUnion::~TGeoVGMultiUnion()
{}

////////////////////////////////////////////////////////////////////////////////
/// compute bounding box of the sphere

void TGeoVGMultiUnion::ComputeBBox()
{
  fDX = fCompositeShape->GetDX();
  fDY = fCompositeShape->GetDY();
  fDZ = fCompositeShape->GetDZ();
}

////////////////////////////////////////////////////////////////////////////////
/// Divide all range of iaxis in range/step cells

TGeoVolume *TGeoVGMultiUnion::Divide(TGeoVolume * /*voldiv*/, const char * /*divname*/, Int_t /*iaxis*/,
                                       Int_t /*ndiv*/, Double_t /*start*/, Double_t /*step*/)
{
   Error("Divide", "Composite shapes cannot be divided");
   return nullptr;
}

////////////////////////////////////////////////////////////////////////////////
/// Returns numbers of vertices, segments and polygons composing the shape mesh.

void TGeoVGMultiUnion::GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const
{
   fCompositeShape->GetMeshNumbers(nvert, nsegs, npols);
}

////////////////////////////////////////////////////////////////////////////////
/// print shape parameters

void TGeoVGMultiUnion::InspectShape() const
{
   printf("*** TGeoVGMultiUnion : %s = %s\n", GetName(), GetTitle());
   printf(" Bounding box:\n");
   fCompositeShape->InspectShape();
}

////////////////////////////////////////////////////////////////////////////////
/// Register the shape and all components to TGeoManager class.

void TGeoVGMultiUnion::RegisterYourself()
{
  fCompositeShape->RegisterYourself();
}

////////////////////////////////////////////////////////////////////////////////
/// Save a primitive as a C++ statement(s) on output stream "out".

void TGeoVGMultiUnion::SavePrimitive(std::ostream &out, Option_t *option /*= ""*/)
{
  fCompositeShape->SavePrimitive(out, option);
}

////////////////////////////////////////////////////////////////////////////////
/// create points for a composite shape

void TGeoVGMultiUnion::SetPoints(Double_t *points) const
{
  fCompositeShape->SetPoints(points);
}

////////////////////////////////////////////////////////////////////////////////
/// create points for a composite shape

void TGeoVGMultiUnion::SetPoints(Float_t *points) const
{
  fCompositeShape->SetPoints(points);
}

////////////////////////////////////////////////////////////////////////////////
/// compute size of this 3D object

void TGeoVGMultiUnion::Sizeof3D() const
{
  fCompositeShape->Sizeof3D();
}

////////////////////////////////////////////////////////////////////////////////
/// Return number of vertices of the mesh representation

Int_t TGeoVGMultiUnion::GetNmeshVertices() const
{
  return fCompositeShape->GetNmeshVertices();
}
