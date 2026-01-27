/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGAdapter
#define ROOT_TGeoVGAdapter

// #if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoBBox.h"

#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>

#include <set>

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoVGAdapter - bridge class for using a VecGeom solid as TGeoShape.   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

namespace vecgeom {
class Transformation3D;
class LogicalVolume;
class VPlacedVolume;
class VUnplacedVolume;
}

template <class UnplacedVolume_t>
class TGeoVGAdapter : public TGeoBBox, protected UnplacedVolume_t {
public:
    using U3Vector = vecgeom::Vector3D<Double_t>;

    /** VecGeom volumes have special delete/new ("AlignedBase")
        and we need to make these functions public again. */
    using UnplacedVolume_t::operator new;
    using UnplacedVolume_t::operator new[];
    using UnplacedVolume_t::operator delete;
    using UnplacedVolume_t::operator delete[];    
    // using TObject::operator delete;
    // using TObject::operator delete[];
    // using TObject::operator new;
    // using TObject::operator new[];

     /**
     * Default Constructor
     */
    // TGeoVGAdapter(const TString& name);

    /**
     * Constructor templated on arguments for UnplacedVolume_t.
     *  @param[in] name The name of the volume.
     *  @param[in] params Templated arguments for UnplacedVolume_t.
     */
    template <typename... T>
    TGeoVGAdapter(const TString& name, const T &... params);

    /**
     * Virtual destructor.
     */
    virtual ~TGeoVGAdapter();

    /**
     * Copy constructor and assignment operator.
     */
    TGeoVGAdapter(const TGeoVGAdapter& rhs) = delete;
    TGeoVGAdapter& operator=(const TGeoVGAdapter& rhs) = delete;

    /**
     * Equality operator. Returns true only if addresses are the same.
     */
    Bool_t operator==(const TGeoVGAdapter& s) const;

    // Overriding functions (from TGeoBBox::TGeoShape)
    //
   // Double_t Capacity() const override; Not needed as done via UnplacedVolume_t
   // void ComputeBBox() override;  Not needed as done via TGeoBBox
   void ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const override;
   Bool_t Contains(const Double_t *point) const override;
   // Bool_t CouldBeCrossed(const Double_t *point, const Double_t *dir) const override
   // {
   //    return fShape->CouldBeCrossed(point, dir);
   // }
   // Int_t DistancetoPrimitive(Int_t px, Int_t py) override { return fShape->DistancetoPrimitive(px, py); }
   Double_t DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact = 1, Double_t step = TGeoShape::Big(),
                           Double_t *safe = nullptr) const override;
   Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                            Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   TGeoVolume *Divide(TGeoVolume *, const char *, Int_t, Int_t, Double_t, Double_t) override { return nullptr; }
   // void Draw(Option_t *option = "") override { fShape->Draw(option); } // *MENU*
   // const char *GetAxisName(Int_t iaxis) const override { return (fShape->GetAxisName(iaxis)); }
   // Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override
   // {
   //    return (fShape->GetAxisRange(iaxis, xlo, xhi));
   // }
   // void GetBoundingCylinder(Double_t *param) const override { return (fShape->GetBoundingCylinder(param)); }
   // const TBuffer3D &GetBuffer3D(Int_t reqSections, Bool_t localFrame) const override
   // {
   //    return (fShape->GetBuffer3D(reqSections, localFrame));
   // }
   // Int_t GetByteCount() const override { return (fShape->GetByteCount()); }
   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override;
   // Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override
   // {
   //    return (fShape->GetPointsOnSegments(npoints, array));
   // }
   // Int_t
   // GetFittingBox(const TGeoBBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const override
   // {
   //    return (fShape->GetFittingBox(parambox, mat, dx, dy, dz));
   // }
   // TGeoShape *GetMakeRuntimeShape(TGeoShape *mother, TGeoMatrix *mat) const override
   // {
   //    return (fShape->GetMakeRuntimeShape(mother, mat));
   // }
   // void GetMeshNumbers(Int_t &nvert, Int_t &nsegs, Int_t &npols) const override
   // {
   //    fShape->GetMeshNumbers(nvert, nsegs, npols);
   // }
   // Not needed as done via TGeoBBox
   // const char *GetName() const override { return (fShape->GetName()); }
   // Int_t GetNmeshVertices() const override { return (fShape->GetNmeshVertices()); }
   // Bool_t IsAssembly() const override { return (fShape->IsAssembly()); }
   // Bool_t IsComposite() const override { return (fShape->IsComposite()); }
   // Bool_t IsCylType() const override { return (fShape->IsCylType()); }
   // Bool_t IsReflected() const override { return (fShape->IsReflected()); }
   // Bool_t IsValidBox() const override { return (fShape->IsValidBox()); }
   // Bool_t IsVecGeom() const override { return kTRUE; }
   void InspectShape() const override;
   // TBuffer3D *MakeBuffer3D() const override { return (fShape->MakeBuffer3D()); }
   // void Paint(Option_t *option = "") override { fShape->Paint(option); }
   // void SetDimensions(Double_t *param) override { fShape->SetDimensions(param); }
   // void SetPoints(Double_t *points) const override { fShape->SetPoints(points); }
   // void SetPoints(Float_t *points) const override { fShape->SetPoints(points); }
   // void SetSegsAndPols(TBuffer3D &buff) const override { fShape->SetSegsAndPols(buff); }
   // void Sizeof3D() const override { fShape->Sizeof3D(); }

   // vector functions
   void Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const override;
   void ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize) override;
   void DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                         Double_t *step) const override;
   void DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                          Double_t *step) const override;
   void Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const override;

    // VecGeom overridden methods ---------------------------------------------

   vecgeom::Precision
   DistanceToOut(U3Vector const& position, U3Vector const& direction,
                 vecgeom::Precision stepMax = vecgeom::kInfinityDbl) const override
   {
     return UnplacedVolume_t::DistanceToOut(position, direction, stepMax);
   }

   vecgeom::EnumInside
   Inside(U3Vector const& aPoint) const override
   {
     return UnplacedVolume_t::Inside(aPoint);
   }

   vecgeom::Precision
   DistanceToIn(U3Vector const& position, U3Vector const& direction,
                const vecgeom::Precision step_max = vecgeom::kInfinityDbl) const override
   {
     return UnplacedVolume_t::DistanceToIn(position, direction, step_max);
   }

   Bool_t Normal(U3Vector const& aPoint, U3Vector& aNormal) const override
   {
     return UnplacedVolume_t::Normal(aPoint, aNormal);
   }

   void Extent(U3Vector& aMin, U3Vector& aMax) const override
   {
     return UnplacedVolume_t::Extent(aMin, aMax);
   }

   U3Vector SamplePointOnSurface() const override
   {
     return UnplacedVolume_t::SamplePointOnSurface();
   }

protected:
   using UnplacedVolume_t::Contains;
   using UnplacedVolume_t::DistanceToOut;
   using UnplacedVolume_t::DistanceToIn;
   using TGeoBBox::Inside;

   ClassDefOverride(TGeoVGAdapter, 1) // Adapter for a VecGeom shape
};

// Inline implementations

#include "TGeoVGAdapter.icc"

// #endif //  ROOT_USE_VECGEOM_SOLIDS

#endif
