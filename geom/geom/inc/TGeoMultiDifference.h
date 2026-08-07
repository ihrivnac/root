// Author: ROOT team

#ifndef ROOT_TGeoMultiDifference
#define ROOT_TGeoMultiDifference

#include "TGeoMultiUnion.h"

#include <vector>

/// A difference between a union of positive nodes and a union of negative nodes.
///
/// All nodes share one acceleration structure. This avoids the duplicate
/// classification, transformations and BVH traversals required by a
/// TGeoSubtraction containing two TGeoMultiUnion shapes.
class TGeoMultiDifference : public TGeoMultiUnion {
private:
   Int_t fNpositive{0};                // Number of positive nodes at the front of the node array
   Bool_t fAddingNegative{kFALSE};     // Whether construction has moved to negative nodes
   std::vector<UChar_t> fBVHSignMask; ///<! Positive/negative membership of every BVH subtree

   using TGeoMultiUnion::AddNode;

   void BuildDifferenceBVH();
   void BuildSignMasks();
   Bool_t GroupContains(const Double_t *point, Bool_t positive) const;
   Double_t GroupDistFromInside(const Double_t *point, const Double_t *dir, Bool_t positive, Double_t step) const;
   Double_t GroupDistFromOutside(const Double_t *point, const Double_t *dir, Bool_t positive, Double_t step) const;

public:
   TGeoMultiDifference();
   explicit TGeoMultiDifference(const char *name);
   ~TGeoMultiDifference() override;

   void AddNegativeNode(TGeoShape &shape, const TGeoMatrix &matrix);
   void AddNegativeNode(TGeoShape *shape, const TGeoMatrix *matrix = nullptr);
   void AddPositiveNode(TGeoShape &shape, const TGeoMatrix &matrix);
   void AddPositiveNode(TGeoShape *shape, const TGeoMatrix *matrix = nullptr);
   void AfterStreamer() override;
   Double_t Capacity() const override;
   void ComputeBBox() override;
   void ComputeNormal(const Double_t *point, const Double_t *dir, Double_t *norm) const override;
   Bool_t Contains(const Double_t *point) const override;
   Double_t DistFromInside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                           Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   Double_t DistFromOutside(const Double_t *point, const Double_t *dir, Int_t iact = 1,
                            Double_t step = TGeoShape::Big(), Double_t *safe = nullptr) const override;
   Int_t GetNnegative() const { return GetNnodes() - fNpositive; }
   Int_t GetNpositive() const { return fNpositive; }
   void InspectShape() const override;
   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override;
   void SavePrimitive(std::ostream &out, Option_t *option = "") override;
   void Voxelize();

   ClassDefOverride(TGeoMultiDifference, 1) // Difference of positive and negative multi-unions
};

#endif
