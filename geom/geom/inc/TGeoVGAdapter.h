/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGAdapter
#define ROOT_TGeoVGAdapter

#if defined(ROOT_USE_VECGEOM_SOLIDS)

#include "TGeoShape.h"
#include "TGeoBaseBox.h"

#include <VecGeom/base/Global.h>
#include <VecGeom/base/Vector3D.h>

#include <iostream>

#ifdef __ROOTCLING__
// Mocking the base class for rootcling to allow inheritance and 'override'
namespace vecgeom {
  class VUnplacedVolume {
  public:
    virtual ~VUnplacedVolume() {}
    
    virtual bool Contains(Vector3D<Precision> const&) const { return false; }
    virtual void Extent(Vector3D<Precision>&, Vector3D<Precision>&) const {}
    virtual bool Normal(Vector3D<Precision> const&, Vector3D<Precision>&) const { return false; }
    virtual Precision DistanceToIn(Vector3D<Precision> const&, Vector3D<Precision> const&, const Precision) const { return 0; }
    virtual Precision DistanceToOut(Vector3D<Precision> const&, Vector3D<Precision> const&, const Precision) const { return 0; }
    virtual Precision SafetyToIn(Vector3D<Precision> const&) const { return 0; }
    virtual Precision SafetyToOut(Vector3D<Precision> const&) const { return 0; }
    virtual EnumInside Inside(Vector3D<Precision> const&) const { return kOutside; }
    virtual Precision Capacity() const { return 0; }
    virtual void Print() const { }
    
    // This was the missing function causing the 'override' error
    virtual Vector3D<Precision> SamplePointOnSurface() const { return Vector3D<Precision>(0,0,0); }
    
    // Satisfy warnings by returning a dummy address instead of nullptr
    static void* operator new(size_t s) noexcept { return (void*)0x11; }
    static void* operator new[](size_t s) noexcept { return (void*)0x11; }
    static void operator delete(void*) noexcept {}
    static void operator delete[](void*) noexcept {}
  };
}
#else
#include <VecGeom/volumes/UnplacedVolume.h>
#endif

////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TGeoVGAdapter - bridge class for using a VecGeom solid as TGeoShape.   //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

// namespace vecgeom {
// class Transformation3D;
// class LogicalVolume;
// class VPlacedVolume;
// class VUnplacedVolume;
// }

template <class UnplacedVolume_t>
class TGeoVGAdapter : public TGeoShape, protected UnplacedVolume_t 
{
protected:
     // Bounding box
     TGeoBaseBox fBoundingBox; // bounding box

public:
    using U3Vector = vecgeom::Vector3D<Double_t>;

    /** VecGeom volumes have special delete/new ("AlignedBase")
        and we need to make these functions public again. */
    using UnplacedVolume_t::operator new;
    using UnplacedVolume_t::operator new[];
    using UnplacedVolume_t::operator delete;
    using UnplacedVolume_t::operator delete[];

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
   Double_t Capacity() const override;
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
   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override;
   void InspectShape() const override;

   // vector functions
   void Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const override;
   void ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize) override;
   void DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                         Double_t *step) const override;
   void DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                          Double_t *step) const override;
   void Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const override;
   Bool_t IsVecGeom() const override { return kTRUE; }

   // Bounding box functions
   const TGeoBaseBox* GetBoundingBox() const override { return &fBoundingBox; }
   virtual Double_t GetDX() const { return fBoundingBox.GetDX(); }
   virtual Double_t GetDY() const { return fBoundingBox.GetDY(); }
   virtual Double_t GetDZ() const { return fBoundingBox.GetDZ(); }
   void SetBoxDimensions(Double_t dx, Double_t dy, Double_t dz, Double_t *origin = nullptr) { fBoundingBox.SetBoxDimensions(dx, dy, dz, origin); }
   void SetDimensions(Double_t *param) override { fBoundingBox.SetDimensions(param); }

   void ComputeBBox() override {fBoundingBox.ComputeBBox(); }
   Bool_t CouldBeCrossed(const Double_t *point, const Double_t *dir) const override { return fBoundingBox.CouldBeCrossed(point, dir); }
   Int_t DistancetoPrimitive(Int_t px, Int_t py) override { return fBoundingBox.DistancetoPrimitive(px, py); }
   const char *GetAxisName(Int_t iaxis) const override  { return fBoundingBox.GetAxisName(iaxis); }
   Double_t GetAxisRange(Int_t iaxis, Double_t &xlo, Double_t &xhi) const override { return fBoundingBox.GetAxisRange(iaxis, xlo, xhi); }
   void GetBoundingCylinder(Double_t *param) const override { fBoundingBox.GetBoundingCylinder(param); }
   Int_t GetByteCount() const override { return fBoundingBox.GetByteCount(); }
   virtual Bool_t GetPointsOnFacet(Int_t index, Int_t npoints, Double_t *array) const { return fBoundingBox.GetPointsOnFacet(index, npoints, array); }
   Bool_t GetPointsOnSegments(Int_t npoints, Double_t *array) const override  { return fBoundingBox.GetPointsOnSegments(npoints, array); }
   Int_t
   GetFittingBox(const TGeoBaseBox *parambox, TGeoMatrix *mat, Double_t &dx, Double_t &dy, Double_t &dz) const override { return fBoundingBox.GetFittingBox(parambox, mat, dx, dy, dz); }
   Bool_t IsCylType() const override { return fBoundingBox.IsCylType(); }
   Bool_t IsValidBox() const override { return fBoundingBox.IsValidBox(); }
   void SetPoints(Double_t *points) const override { return fBoundingBox.SetPoints(points); }
   void SetPoints(Float_t *points) const override { return fBoundingBox.SetPoints(points); }
   void SetSegsAndPols(TBuffer3D &buffer) const override { fBoundingBox.SetSegsAndPols(buffer); }
   void Sizeof3D() const override { fBoundingBox.Sizeof3D(); }

    // VecGeom overridden methods ---------------------------------------------

   vecgeom::Precision
   DistanceToOut(U3Vector const& position, U3Vector const& direction,
                 vecgeom::Precision stepMax = vecgeom::kInfinityDbl) const override
   {
     std::cout << "Calling UnplacedVolume_t DistanceToOut" << std::endl;
     return UnplacedVolume_t::DistanceToOut(position, direction, stepMax);
   }

   vecgeom::EnumInside
   Inside(U3Vector const& aPoint) const override
   {
     std::cout << "Calling UnplacedVolume_t Inside" << std::endl;
     return UnplacedVolume_t::Inside(aPoint);
   }

   vecgeom::Precision
   DistanceToIn(U3Vector const& position, U3Vector const& direction,
                const vecgeom::Precision step_max = vecgeom::kInfinityDbl) const override
   {
     std::cout << "Calling UnplacedVolume_t DistanceToIn" << std::endl;
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
   using UnplacedVolume_t::Capacity;
   using UnplacedVolume_t::Contains;
   using UnplacedVolume_t::DistanceToOut;
   using UnplacedVolume_t::DistanceToIn;
   using UnplacedVolume_t::Inside;
   using TGeoShape::Inside;

   ClassDefOverride(TGeoVGAdapter, 1) // Adapter for a VecGeom shape
};

// Inline implementations

#include "TGeoVGAdapter.icc"

#endif //  ROOT_USE_VECGEOM_SOLIDS

#endif
