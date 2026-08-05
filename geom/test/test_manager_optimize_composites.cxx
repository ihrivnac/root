#include <gtest/gtest.h>

#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoMultiUnion.h>
#include <TGeoVolume.h>
#include <TGeoVoxelFinder.h>

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

   auto *ineligibleShape =
      new TGeoCompositeShape("placed_intersection", new TGeoIntersection(first, second));
   auto *ineligibleVolume = new TGeoVolume("ineligible_volume", ineligibleShape, medium);
   top->AddNode(ineligibleVolume, 1);

   auto *unusedUnion2 = new TGeoCompositeShape("unused_union2", new TGeoUnion(first, second));
   auto *unusedShape = new TGeoCompositeShape("unused_union3", new TGeoUnion(unusedUnion2, third));
   auto *unusedVolume = new TGeoVolume("unused_volume", unusedShape, medium);

   manager.CloseGeometry();
   const Int_t shapesBefore = manager.GetListOfShapes()->GetEntriesFast();

   EXPECT_EQ(manager.OptimizeCompositeShapes(), 1);
   EXPECT_EQ(manager.GetListOfShapes()->GetEntriesFast(), shapesBefore);
   EXPECT_EQ(firstVolume->GetShape(), placedShape);
   EXPECT_EQ(secondVolume->GetShape(), placedShape);
   EXPECT_EQ(unusedVolume->GetShape(), unusedShape);

   ASSERT_NE(top->GetVoxels(), nullptr);
   EXPECT_FALSE(top->GetVoxels()->NeedRebuild());
   EXPECT_EQ(manager.OptimizeCompositeShapes(kTRUE), 1);

   auto *replacement = dynamic_cast<TGeoMultiUnion *>(firstVolume->GetShape());
   ASSERT_NE(replacement, nullptr);
   EXPECT_EQ(replacement->GetNumberOfSolids(), 4);
   EXPECT_EQ(secondVolume->GetShape(), replacement);
   EXPECT_EQ(ineligibleVolume->GetShape(), ineligibleShape);
   EXPECT_EQ(unusedVolume->GetShape(), unusedShape);
   EXPECT_TRUE(top->GetVoxels()->NeedRebuild());
}
