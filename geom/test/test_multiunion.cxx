#include <gtest/gtest.h>

#include <TBufferFile.h>
#include <TClass.h>
#include <TGeoBBox.h>
#include <TGeoMatrix.h>
#include <TGeoMultiUnion.h>

#include <memory>

namespace {
constexpr double kTol = 1.e-9;
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

TEST(TGeoMultiUnion, StreamingRebuildsAccelerationData)
{
   TGeoBBox box(1., 1., 1.);
   TGeoTranslation transform(3., 0., 0.);
   TGeoMultiUnion original;
   original.AddNode(box, transform);
   original.Voxelize();

   TBufferFile buffer(TBuffer::kWrite);
   buffer.WriteObjectAny(&original, TGeoMultiUnion::Class());
   buffer.SetReadMode();
   buffer.SetBufferOffset(0);
   std::unique_ptr<TGeoMultiUnion> restored(
      static_cast<TGeoMultiUnion *>(buffer.ReadObjectAny(TGeoMultiUnion::Class())));

   ASSERT_NE(restored, nullptr);
   restored->AfterStreamer();
   EXPECT_TRUE(restored->IsVoxelized());
   EXPECT_EQ(restored->GetNnodes(), 1);
   const double center[3] = {3., 0., 0.};
   EXPECT_TRUE(restored->Contains(center));
}
