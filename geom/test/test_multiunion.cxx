#include <gtest/gtest.h>

#include <TBufferFile.h>
#include <TClass.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoMultiUnion.h>

#include <memory>
#include <vector>

namespace {
constexpr double kTol = 1.e-9;

class CountingBox : public TGeoBBox {
public:
   using TGeoBBox::TGeoBBox;

   Bool_t Contains(const Double_t *point) const override
   {
      ++fContainsCalls;
      return TGeoBBox::Contains(point);
   }

   Double_t Safety(const Double_t *point, Bool_t in = kTRUE) const override
   {
      ++fSafetyCalls;
      return TGeoBBox::Safety(point, in);
   }

   void ResetCounts() const
   {
      fContainsCalls = 0;
      fSafetyCalls = 0;
   }

   mutable int fContainsCalls = 0;
   mutable int fSafetyCalls = 0;
};
} // namespace

TEST(TGeoMultiUnion, DisjointNodes)
{
   TGeoBBox left(1., 2., 3.);
   TGeoBBox right(1., 2., 3.);
   TGeoTranslation leftTransform(-4., 0., 0.);
   TGeoTranslation rightTransform(4., 0., 0.);
   TGeoMultiUnion shape("disjoint");
   shape.AddNode(left, leftTransform);
   shape.AddNode(right, rightTransform);
   shape.Voxelize();

   EXPECT_TRUE(shape.IsVoxelized());
   EXPECT_EQ(shape.GetNumberOfSolids(), 2);
   const double inLeft[3] = {-4., 0., 0.};
   const double inRight[3] = {4., 0., 0.};
   const double between[3] = {0., 0., 0.};
   EXPECT_TRUE(shape.Contains(inLeft));
   EXPECT_TRUE(shape.Contains(inRight));
   EXPECT_FALSE(shape.Contains(between));

   const double fromOutside[3] = {-10., 0., 0.};
   const double plusX[3] = {1., 0., 0.};
   EXPECT_NEAR(shape.DistFromOutside(fromOutside, plusX, 3), 5., kTol);
   EXPECT_NEAR(shape.DistFromInside(inLeft, plusX, 3), 1., kTol);

   EXPECT_NEAR(shape.GetDX(), 5., kTol);
   EXPECT_NEAR(shape.GetDY(), 2., kTol);
   EXPECT_NEAR(shape.GetDZ(), 3., kTol);
}

TEST(TGeoMultiUnion, OverlapIsNavigatedAsOneSolid)
{
   TGeoBBox first(2., 1., 1.);
   TGeoBBox second(2., 1., 1.);
   TGeoTranslation firstTransform(-1., 0., 0.);
   TGeoTranslation secondTransform(1., 0., 0.);
   TGeoMultiUnion shape;
   shape.AddNode(first, firstTransform);  // [-3, 1]
   shape.AddNode(second, secondTransform); // [-1, 3]
   shape.Voxelize();

   const double center[3] = {0., 0., 0.};
   const double plusX[3] = {1., 0., 0.};
   const double minusX[3] = {-1., 0., 0.};
   EXPECT_NEAR(shape.DistFromInside(center, plusX, 3), 3., kTol);
   EXPECT_NEAR(shape.DistFromInside(center, minusX, 3), 3., kTol);

   const double leftCenter[3] = {-2., 0., 0.};
   EXPECT_NEAR(shape.DistFromInside(leftCenter, plusX, 3), 5., 5.e-9);
}

TEST(TGeoMultiUnion, TouchingNodesHaveContinuousRayNavigation)
{
   TGeoBBox first(1., 1., 1.);
   TGeoBBox second(1., 1., 1.);
   TGeoTranslation firstTransform(-1., 0., 0.);
   TGeoTranslation secondTransform(1., 0., 0.);
   TGeoMultiUnion shape;
   shape.AddNode(first, firstTransform);   // [-2, 0]
   shape.AddNode(second, secondTransform); // [0, 2]
   shape.Voxelize();

   const double point[3] = {-1., 0., 0.};
   const double direction[3] = {1., 0., 0.};
   EXPECT_NEAR(shape.DistFromInside(point, direction, 3), 3., 5.e-9);
}

TEST(TGeoMultiUnion, RotatedNode)
{
   TGeoBBox box(2., 1., 1.);
   TGeoRotation rotation;
   rotation.RotateZ(90.);
   TGeoCombiTrans transform(3., 0., 0., &rotation);
   TGeoMultiUnion shape;
   shape.AddNode(box, transform);
   shape.Voxelize();

   const double alongRotatedLongAxis[3] = {3., 1.5, 0.};
   const double outsideShortAxis[3] = {4.5, 0., 0.};
   EXPECT_TRUE(shape.Contains(alongRotatedLongAxis));
   EXPECT_FALSE(shape.Contains(outsideShortAxis));
}

TEST(TGeoMultiUnion, SafetyClassifiesEachNodeOnlyOnce)
{
   CountingBox left(1., 1., 1.);
   CountingBox right(1., 1., 1.);
   TGeoTranslation leftTransform(-4., 0., 0.);
   TGeoTranslation rightTransform(4., 0., 0.);
   TGeoMultiUnion shape;
   shape.AddNode(left, leftTransform);
   shape.AddNode(right, rightTransform);
   shape.Voxelize();

   const double insideLeft[3] = {-4., 0., 0.};
   EXPECT_NEAR(shape.Safety(insideLeft, kTRUE), 1., kTol);
   EXPECT_EQ(left.fContainsCalls, 1);
   EXPECT_EQ(left.fSafetyCalls, 1);
   EXPECT_EQ(right.fContainsCalls, 0);
   EXPECT_EQ(right.fSafetyCalls, 0);

   left.ResetCounts();
   right.ResetCounts();
   const double outside[3] = {0., 0., 0.};
   EXPECT_NEAR(shape.Safety(outside, kFALSE), 3., kTol);
   EXPECT_EQ(left.fContainsCalls, 0);
   EXPECT_EQ(left.fSafetyCalls, 1);
   EXPECT_EQ(right.fContainsCalls, 0);
   EXPECT_EQ(right.fSafetyCalls, 1);
}

TEST(TGeoMultiUnion, SmallUnionSafetyUsesNodeBoxLowerBounds)
{
   CountingBox near(1., 1., 1.);
   CountingBox far(1., 1., 1.);
   TGeoTranslation nearTransform(4., 0., 0.);
   TGeoTranslation farTransform(100., 0., 0.);
   TGeoMultiUnion shape;
   shape.AddNode(near, nearTransform);
   shape.AddNode(far, farTransform);
   shape.Voxelize();

   const double outside[3] = {0., 0., 0.};
   EXPECT_NEAR(shape.Safety(outside, kFALSE), 3., kTol);
   EXPECT_EQ(near.fContainsCalls, 0);
   EXPECT_EQ(near.fSafetyCalls, 1);
   EXPECT_EQ(far.fContainsCalls, 0);
   EXPECT_EQ(far.fSafetyCalls, 0);
}

TEST(TGeoMultiUnion, LargeUnionUsesBVH)
{
   std::vector<std::unique_ptr<CountingBox>> boxes;
   std::vector<std::unique_ptr<TGeoTranslation>> transforms;
   TGeoMultiUnion shape;
   for (int index = 0; index < 16; ++index) {
      boxes.emplace_back(std::make_unique<CountingBox>(1., 1., 1.));
      transforms.emplace_back(std::make_unique<TGeoTranslation>(10. * (index + 1), 0., 0.));
      shape.AddNode(*boxes.back(), *transforms.back());
   }
   EXPECT_FALSE(shape.IsBVHEnabled());
   shape.Voxelize();

   EXPECT_TRUE(shape.IsBVHEnabled());
   const double outside[3] = {0., 0., 0.};
   EXPECT_NEAR(shape.Safety(outside, kFALSE), 9., kTol);
   int safetyCalls = 0;
   int containsCalls = 0;
   for (const auto &box : boxes) {
      safetyCalls += box->fSafetyCalls;
      containsCalls += box->fContainsCalls;
   }
   EXPECT_EQ(safetyCalls, 1);
   EXPECT_EQ(containsCalls, 0);

   for (const auto &box : boxes)
      box->ResetCounts();
   const double insideLast[3] = {160., 0., 0.};
   EXPECT_TRUE(shape.Contains(insideLast));
   containsCalls = 0;
   for (const auto &box : boxes)
      containsCalls += box->fContainsCalls;
   EXPECT_EQ(containsCalls, 1);

   const double plusX[3] = {1., 0., 0.};
   EXPECT_NEAR(shape.DistFromOutside(outside, plusX, 3), 9., kTol);
}

TEST(TGeoMultiUnion, SafetyRejectsIncorrectInsideStateInOnePass)
{
   CountingBox first(1., 1., 1.);
   CountingBox second(1., 1., 1.);
   TGeoTranslation firstTransform(-4., 0., 0.);
   TGeoTranslation secondTransform(4., 0., 0.);
   TGeoMultiUnion shape;
   shape.AddNode(first, firstTransform);
   shape.AddNode(second, secondTransform);

   const double insideFirst[3] = {-4., 0., 0.};
   EXPECT_DOUBLE_EQ(shape.Safety(insideFirst, kFALSE), 0.);
   EXPECT_EQ(first.fContainsCalls, 1);
   EXPECT_EQ(first.fSafetyCalls, 0);
   EXPECT_EQ(second.fContainsCalls, 0);
   EXPECT_EQ(second.fSafetyCalls, 0);

   first.ResetCounts();
   second.ResetCounts();
   const double outside[3] = {0., 0., 0.};
   EXPECT_DOUBLE_EQ(shape.Safety(outside, kTRUE), 0.);
   EXPECT_EQ(first.fContainsCalls, 1);
   EXPECT_EQ(first.fSafetyCalls, 0);
   EXPECT_EQ(second.fContainsCalls, 1);
   EXPECT_EQ(second.fSafetyCalls, 0);
}

TEST(TGeoMultiUnion, StreamingRebuildsAccelerationData)
{
   std::vector<std::unique_ptr<TGeoBBox>> boxes;
   std::vector<std::unique_ptr<TGeoTranslation>> transforms;
   TGeoMultiUnion original;
   for (int index = 0; index < 12; ++index) {
      boxes.emplace_back(std::make_unique<TGeoBBox>(1., 1., 1.));
      transforms.emplace_back(std::make_unique<TGeoTranslation>(3. * index, 0., 0.));
      original.AddNode(*boxes.back(), *transforms.back());
   }
   original.Voxelize();
   ASSERT_TRUE(original.IsBVHEnabled());

   TBufferFile buffer(TBuffer::kWrite);
   buffer.WriteObjectAny(&original, TGeoMultiUnion::Class());
   buffer.SetReadMode();
   buffer.SetBufferOffset(0);
   std::unique_ptr<TGeoMultiUnion> restored(
      static_cast<TGeoMultiUnion *>(buffer.ReadObjectAny(TGeoMultiUnion::Class())));

   ASSERT_NE(restored, nullptr);
   restored->AfterStreamer();
   EXPECT_TRUE(restored->IsVoxelized());
   EXPECT_TRUE(restored->IsBVHEnabled());
   EXPECT_EQ(restored->GetNnodes(), 12);
   const double center[3] = {33., 0., 0.};
   EXPECT_TRUE(restored->Contains(center));
}
