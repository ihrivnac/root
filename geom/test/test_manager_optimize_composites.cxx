#include <gtest/gtest.h>

#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoMultiDifference.h>
#include <TGeoMultiUnion.h>
#include <TGeoVolume.h>
#include <TGeoVoxelFinder.h>

#include <string>

namespace {

TGeoCompositeShape *MakeNestedUnion(const char *prefix, Int_t leaves)
{
   const auto makeBox = [prefix](Int_t index) {
      const std::string name = std::string(prefix) + "_box_" + std::to_string(index);
      return new TGeoBBox(name.c_str(), 1., 1., 1.);
   };
   std::string name = std::string(prefix) + "_union_2";
   auto *result = new TGeoCompositeShape(name.c_str(), new TGeoUnion(makeBox(0), makeBox(1)));
   for (Int_t index = 2; index < leaves; ++index) {
      name = std::string(prefix) + "_union_" + std::to_string(index + 1);
      result = new TGeoCompositeShape(name.c_str(), new TGeoUnion(result, makeBox(index)));
   }
   return result;
}

} // namespace

TEST(TGeoManagerOptimizeCompositeShapes, ProcessesOnlyDistinctPlacedShapes)
{
   TGeoManager manager("manager_optimize", "manager optimize test");
   auto *material = new TGeoMaterial("material", 1., 1., 1.);
   auto *medium = new TGeoMedium("medium", 1, material);
   auto *top = new TGeoVolume("top", new TGeoBBox("top_box", 100., 100., 100.), medium);
   manager.SetTopVolume(top);

   auto *first = new TGeoBBox("first", 1., 1., 1.);
   auto *second = new TGeoBBox("second", 1., 1., 1.);
   auto *third = new TGeoBBox("third", 1., 1., 1.);
   auto *fourth = new TGeoBBox("fourth", 1., 1., 1.);
   auto *firstMatrix = new TGeoTranslation("first_matrix", -3., 0., 0.);
   auto *secondMatrix = new TGeoTranslation("second_matrix", -1., 0., 0.);
   auto *thirdMatrix = new TGeoTranslation("third_matrix", 1., 0., 0.);
   auto *fourthMatrix = new TGeoTranslation("fourth_matrix", 3., 0., 0.);

   auto *union2 = new TGeoCompositeShape("component_union2",
                                         new TGeoUnion(first, second, firstMatrix, secondMatrix));
   auto *componentOnly = new TGeoCompositeShape("component_union3",
                                                 new TGeoUnion(union2, third, nullptr, thirdMatrix));
   auto *placedShape = new TGeoCompositeShape("placed_union4",
                                               new TGeoUnion(componentOnly, fourth, nullptr, fourthMatrix));

   auto *firstVolume = new TGeoVolume("first_volume", placedShape, medium);
   auto *secondVolume = new TGeoVolume("second_volume", placedShape, medium);
   top->AddNode(firstVolume, 1, new TGeoTranslation(-20., 0., 0.));
   top->AddNode(secondVolume, 1, new TGeoTranslation(20., 0., 0.));

   auto *intersectionFirst = new TGeoBBox("intersection_first", 1., 1., 1.);
   auto *intersectionSecond = new TGeoBBox("intersection_second", 1., 1., 1.);
   auto *ineligibleShape = new TGeoCompositeShape("placed_intersection",
                                                   new TGeoIntersection(intersectionFirst, intersectionSecond));
   auto *ineligibleVolume = new TGeoVolume("ineligible_volume", ineligibleShape, medium);
   top->AddNode(ineligibleVolume, 1);

   auto *unusedUnion2 = new TGeoCompositeShape("unused_union2", new TGeoUnion(first, second));
   auto *unusedShape = new TGeoCompositeShape("unused_union3", new TGeoUnion(unusedUnion2, third));
   auto *unusedVolume = new TGeoVolume("unused_volume", unusedShape, medium);

   manager.CloseGeometry();
   const Int_t shapesBefore = manager.GetListOfShapes()->GetEntriesFast();

   EXPECT_EQ(manager.OptimizeCompositeShapes(kFALSE, 3), 1);
   EXPECT_EQ(manager.GetListOfShapes()->GetEntriesFast(), shapesBefore);
   EXPECT_EQ(firstVolume->GetShape(), placedShape);
   EXPECT_EQ(secondVolume->GetShape(), placedShape);
   EXPECT_EQ(unusedVolume->GetShape(), unusedShape);

   ASSERT_NE(top->GetVoxels(), nullptr);
   EXPECT_FALSE(top->GetVoxels()->NeedRebuild());
   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE, 3), 1);

   auto *replacement = dynamic_cast<TGeoMultiUnion *>(firstVolume->GetShape());
   ASSERT_NE(replacement, nullptr);
   EXPECT_EQ(replacement->GetNumberOfSolids(), 4);
   EXPECT_EQ(secondVolume->GetShape(), replacement);
   EXPECT_EQ(ineligibleVolume->GetShape(), ineligibleShape);
   EXPECT_EQ(unusedVolume->GetShape(), unusedShape);
   EXPECT_TRUE(top->GetVoxels()->NeedRebuild());
}

TEST(TGeoManagerOptimizeCompositeShapes, AppliesAutomaticMinimumLeafThreshold)
{
   TGeoManager manager("manager_threshold", "manager automatic threshold test");
   auto *material = new TGeoMaterial("threshold_material", 1., 1., 1.);
   auto *medium = new TGeoMedium("threshold_medium", 1, material);
   auto *top = new TGeoVolume("threshold_top", new TGeoBBox("threshold_top_box", 100., 100., 100.), medium);
   manager.SetTopVolume(top);

   auto *shallow = MakeNestedUnion("shallow", 4);
   auto *deep = MakeNestedUnion("deep", 10);
   auto *shallowVolume = new TGeoVolume("shallow_volume", shallow, medium);
   auto *deepVolume = new TGeoVolume("deep_volume", deep, medium);
   top->AddNode(shallowVolume, 1, new TGeoTranslation(-10., 0., 0.));
   top->AddNode(deepVolume, 1, new TGeoTranslation(10., 0., 0.));
   manager.CloseGeometry();

   EXPECT_EQ(manager.OptimizeCompositeShapes(), 1);
   EXPECT_EQ(shallowVolume->GetShape(), shallow);
   EXPECT_EQ(deepVolume->GetShape(), deep);

   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE), 1);
   EXPECT_EQ(shallowVolume->GetShape(), shallow);
   auto *deepReplacement = dynamic_cast<TGeoMultiUnion *>(deepVolume->GetShape());
   ASSERT_NE(deepReplacement, nullptr);
   EXPECT_EQ(deepReplacement->GetNumberOfSolids(), 10);

   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE, 3), 1);
   auto *shallowReplacement = dynamic_cast<TGeoMultiUnion *>(shallowVolume->GetShape());
   ASSERT_NE(shallowReplacement, nullptr);
   EXPECT_EQ(shallowReplacement->GetNumberOfSolids(), 4);
}

TEST(TGeoManagerOptimizeCompositeShapes, AppliesIndependentMultiDifferenceThreshold)
{
   TGeoManager manager("manager_mixed_threshold", "manager mixed threshold test");
   auto *material = new TGeoMaterial("mixed_threshold_material", 1., 1., 1.);
   auto *medium = new TGeoMedium("mixed_threshold_medium", 1, material);
   auto *top =
      new TGeoVolume("mixed_threshold_top", new TGeoBBox("mixed_threshold_top_box", 100., 100., 100.), medium);
   manager.SetTopVolume(top);

   auto *ordinaryUnion = MakeNestedUnion("small_ordinary_union", 3);

   auto *positiveUnion = MakeNestedUnion("small_positive_union", 2);
   auto *negativeLeaf = new TGeoBBox("small_negative_leaf", 1., 1., 1.);
   auto *multiDifference =
      new TGeoCompositeShape("small_multi_difference", new TGeoSubtraction(positiveUnion, negativeLeaf));

   auto *positiveLeaf = new TGeoBBox("small_positive_leaf", 1., 1., 1.);
   auto *negativeUnion = MakeNestedUnion("small_negative_union", 2);
   auto *singleMinusUnion =
      new TGeoCompositeShape("small_single_minus_union", new TGeoSubtraction(positiveLeaf, negativeUnion));

   auto *ordinaryVolume = new TGeoVolume("small_ordinary_volume", ordinaryUnion, medium);
   auto *multiDifferenceVolume = new TGeoVolume("small_multi_difference_volume", multiDifference, medium);
   auto *singleMinusUnionVolume = new TGeoVolume("small_single_minus_union_volume", singleMinusUnion, medium);
   top->AddNode(ordinaryVolume, 1, new TGeoTranslation(-20., 0., 0.));
   top->AddNode(multiDifferenceVolume, 1);
   top->AddNode(singleMinusUnionVolume, 1, new TGeoTranslation(20., 0., 0.));
   manager.CloseGeometry();

   EXPECT_EQ(manager.OptimizeCompositeShapes(kFALSE, 6), 0);
   EXPECT_EQ(manager.OptimizeCompositeShapes(kFALSE, 6, 3), 1);
   EXPECT_EQ(ordinaryVolume->GetShape(), ordinaryUnion);
   EXPECT_EQ(multiDifferenceVolume->GetShape(), multiDifference);
   EXPECT_EQ(singleMinusUnionVolume->GetShape(), singleMinusUnion);

   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE, 6, 3), 1);
   auto *replacement = dynamic_cast<TGeoMultiDifference *>(multiDifferenceVolume->GetShape());
   ASSERT_NE(replacement, nullptr);
   EXPECT_EQ(replacement->GetNpositive(), 2);
   EXPECT_EQ(replacement->GetNnegative(), 1);
   EXPECT_EQ(ordinaryVolume->GetShape(), ordinaryUnion);
   EXPECT_EQ(singleMinusUnionVolume->GetShape(), singleMinusUnion);

   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE, 6, 3), 0);
}
