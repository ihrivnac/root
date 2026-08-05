// Author: ROOT team

#ifndef ROOT_TGeoMultiUnion
#define ROOT_TGeoMultiUnion

#include "TGeoBBox.h"
#include "TObjArray.h"

#include <vector>

class TGeoMatrix;

/// A Boolean union of an arbitrary number of transformed TGeoShape objects.
///
/// The interface follows G4MultiUnion: nodes are added one at a time and
/// Voxelize() closes the construction phase and builds per-node bounding boxes
/// used to reject irrelevant shapes during navigation. Constituent shapes are
/// not owned. Transformation matrices are copied and owned by the multi-union.
class TGeoMultiUnion : public TGeoBBox {
private:
   TObjArray fShapes;              // Constituent shapes (not owned)
   TObjArray fMatrices;            //-> Constituent transformations (owned)
   Bool_t fVoxelized{kFALSE};      // Whether the acceleration data is current
   std::vector<Double_t> fBoxes;   ///<! Node AABBs: xmin,xmax,ymin,ymax,zmin,zmax

   TGeoMultiUnion(const TGeoMultiUnion &) = delete;
   TGeoMultiUnion &operator=(const TGeoMultiUnion &) = delete;

   Bool_t AcceptNode(Int_t inode, const Double_t *point) const;
   Bool_t CrossesNodeBox(Int_t inode, const Double_t *point, const Double_t *dir, Double_t step) const;
   void TransformToNode(Int_t inode, const Double_t *point, const Double_t *dir, Double_t *local,
                        Double_t *localDir) const;

public:
   TGeoMultiUnion();
   explicit TGeoMultiUnion(const char *name);
   ~TGeoMultiUnion() override;

   void AddNode(TGeoShape &shape, const TGeoMatrix &matrix);
   void AddNode(TGeoShape *shape, const TGeoMatrix *matrix = nullptr);
   void AfterStreamer() override;
   Double_t Capacity() const override;
   void ClearThreadData() const override;
   void ComputeBBox() override;
   void ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const override;
   void ComputeNormal_v(const Double_t *points, const Double_t *dirs, Double_t *norms, Int_t vecsize) override;
   Bool_t Contains(const Double_t *point) const override;
   void Contains_v(const Double_t *points, Bool_t *inside, Int_t vecsize) const override;
   void CreateThreadData(Int_t nthreads) override;
   Int_t DistancetoPrimitive(Int_t px, Int_t py) override;
   Double_t DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   void DistFromInside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                         Double_t *step) const override;
   Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                            Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   void DistFromOutside_v(const Double_t *points, const Double_t *dirs, Double_t *dists, Int_t vecsize,
                          Double_t *step) const override;
   TGeoVolume *Divide(TGeoVolume *voldiv, const char *divname, Int_t iaxis, Int_t ndiv, Double_t start,
                      Double_t step) override;
   TGeoMatrix *GetMatrix(Int_t inode) const;
   Int_t GetNnodes() const { return fShapes.GetEntriesFast(); }
   Int_t GetNumberOfSolids() const { return GetNnodes(); }
   TGeoShape *GetShape(Int_t inode) const;
   TGeoShape *GetSolid(Int_t inode) const { return GetShape(inode); }
   TGeoMatrix *GetTransformation(Int_t inode) const { return GetMatrix(inode); }
   TGeoShape *GetMakeRuntimeShape(TGeoShape *, TGeoMatrix *) const override { return nullptr; }
   void InspectShape() const override;
   Bool_t IsConvex() const final { return kFALSE; }
   Bool_t IsCylType() const override { return kFALSE; }
   Bool_t IsVoxelized() const { return fVoxelized; }
   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override;
   void Safety_v(const Double_t *points, const Bool_t *inside, Double_t *safe, Int_t vecsize) const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void SetDimensions(Double_t *) override {}
   void Voxelize();

   ClassDefOverride(TGeoMultiUnion, 1) // Multi-union shape
};

#endif
