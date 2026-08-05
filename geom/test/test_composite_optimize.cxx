#include <gtest/gtest.h>

#include <TGeoBBox.h>
#include <TGeoBoolNode.h>
#include <TGeoCompositeShape.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoMultiUnion.h>

namespace {

void ExpectSameContainment(const TGeoShape &first, const TGeoShape &second)
{
   for (double x = -8.; x <= 8.; x += .5) {
      for (double y = -3.; y <= 3.; y += .5) {
         const double point[3] = {x, y, .25};
         EXPECT_EQ(first.Contains(point), second.Contains(point)) << "point = (" << x << ", " << y << ", .25)";
      }
   }
}

} // namespace

TEST(TGeoCompositeShapeOptimize, LeavesTwoComponentShapeUnchanged)
{
   TGeoBBox left(1., 1., 1.);
   TGeoBBox right(1., 1., 1.);
   TGeoCompositeShape shape("two", new TGeoUnion(&left, &right));

   EXPECT_EQ(shape.Optimize(), &shape);
}

TEST(TGeoCompositeShapeOptimize, FlattensUnionAndComposesTransforms)
{
   TGeoBBox first(1., 1., 1.);
   TGeoBBox second(1., 1., 1.);
   TGeoBBox third(1., 1., 1.);
   auto *firstMatrix = new TGeoTranslation(-2., 0., 0.);
   auto *secondMatrix = new TGeoTranslation(2., 0., 0.);
   auto *nestedMatrix = new TGeoTranslation(3., 0., 0.);
   auto *thirdMatrix = new TGeoTranslation(-5., 0., 0.);
   TGeoCompositeShape nested("nested_union", new TGeoUnion(&first, &second, firstMatrix, secondMatrix));
   TGeoCompositeShape original("union3", new TGeoUnion(&nested, &third, nestedMatrix, thirdMatrix));

   TGeoShape *optimized = original.Optimize();
   auto *multiUnion = dynamic_cast<TGeoMultiUnion *>(optimized);
   ASSERT_NE(multiUnion, nullptr);
   EXPECT_EQ(multiUnion->GetNumberOfSolids(), 3);
   for (int i = 0; i < multiUnion->GetNumberOfSolids(); ++i)
      EXPECT_FALSE(multiUnion->GetSolid(i)->IsComposite());
   ExpectSameContainment(original, *optimized);
}

TEST(TGeoCompositeShapeOptimize, MakesShapeMinusMultiUnion)
{
   TGeoBBox outer(6., 2., 2.);
   TGeoBBox firstHole(1., 1., 3.);
   TGeoBBox secondHole(1., 1., 3.);
   auto *outerMatrix = new TGeoTranslation(1., 0., 0.);
   auto *firstMatrix = new TGeoTranslation(-2., 0., 0.);
   auto *secondMatrix = new TGeoTranslation(2., 0., 0.);
   TGeoCompositeShape holes("holes", new TGeoUnion(&firstHole, &secondHole, firstMatrix, secondMatrix));
   TGeoCompositeShape original("shape_minus_union", new TGeoSubtraction(&outer, &holes, outerMatrix));

   TGeoShape *optimized = original.Optimize();
   auto *subtraction = dynamic_cast<TGeoCompositeShape *>(optimized);
   ASSERT_NE(subtraction, nullptr);
   ASSERT_NE(subtraction, &original);
   ASSERT_EQ(subtraction->GetBoolNode()->GetBooleanOperator(), TGeoBoolNode::kGeoSubtraction);
   EXPECT_EQ(subtraction->GetBoolNode()->GetLeftShape(), &outer);
   ASSERT_NE(subtraction->GetBoolNode()->GetLeftMatrix(), nullptr);
   EXPECT_DOUBLE_EQ(subtraction->GetBoolNode()->GetLeftMatrix()->GetTranslation()[0], 1.);
   auto *negative = dynamic_cast<TGeoMultiUnion *>(subtraction->GetBoolNode()->GetRightShape());
   ASSERT_NE(negative, nullptr);
   EXPECT_EQ(negative->GetNumberOfSolids(), 2);
   ExpectSameContainment(original, *optimized);
}

TEST(TGeoCompositeShapeOptimize, MakesDifferenceOfMultiUnions)
{
   TGeoBBox first(2., 2., 2.);
   TGeoBBox second(2., 2., 2.);
   TGeoBBox firstHole(.5, 1., 3.);
   TGeoBBox secondHole(.5, 1., 3.);
   auto *leftMatrix = new TGeoTranslation(-3., 0., 0.);
   auto *rightMatrix = new TGeoTranslation(3., 0., 0.);
   TGeoCompositeShape positive("positive_union", new TGeoUnion(&first, &second, leftMatrix, rightMatrix));
   TGeoCompositeShape negative("negative_union", new TGeoUnion(&firstHole, &secondHole, leftMatrix, rightMatrix));
   TGeoCompositeShape original("union_minus_union", new TGeoSubtraction(&positive, &negative));

   TGeoShape *optimized = original.Optimize();
   auto *subtraction = dynamic_cast<TGeoCompositeShape *>(optimized);
   ASSERT_NE(subtraction, nullptr);
   ASSERT_NE(subtraction, &original);
   auto *positiveMultiUnion = dynamic_cast<TGeoMultiUnion *>(subtraction->GetBoolNode()->GetLeftShape());
   auto *negativeMultiUnion = dynamic_cast<TGeoMultiUnion *>(subtraction->GetBoolNode()->GetRightShape());
   ASSERT_NE(positiveMultiUnion, nullptr);
   ASSERT_NE(negativeMultiUnion, nullptr);
   EXPECT_EQ(positiveMultiUnion->GetNumberOfSolids(), 2);
   EXPECT_EQ(negativeMultiUnion->GetNumberOfSolids(), 2);
   ExpectSameContainment(original, *optimized);
}

TEST(TGeoCompositeShapeOptimize, FlattensChainedSubtractions)
{
   TGeoBBox outer(6., 2., 2.);
   TGeoBBox firstHole(1., 1., 3.);
   TGeoBBox secondHole(1., 1., 3.);
   auto *firstMatrix = new TGeoTranslation(-2., 0., 0.);
   auto *secondMatrix = new TGeoTranslation(2., 0., 0.);
   TGeoCompositeShape firstSubtraction("first_subtraction",
                                       new TGeoSubtraction(&outer, &firstHole, nullptr, firstMatrix));
   TGeoCompositeShape original("chained_subtraction",
                               new TGeoSubtraction(&firstSubtraction, &secondHole, nullptr, secondMatrix));

   TGeoShape *optimized = original.Optimize();
   auto *subtraction = dynamic_cast<TGeoCompositeShape *>(optimized);
   ASSERT_NE(subtraction, nullptr);
   ASSERT_NE(subtraction, &original);
   EXPECT_EQ(subtraction->GetBoolNode()->GetLeftShape(), &outer);
   auto *negative = dynamic_cast<TGeoMultiUnion *>(subtraction->GetBoolNode()->GetRightShape());
   ASSERT_NE(negative, nullptr);
   EXPECT_EQ(negative->GetNumberOfSolids(), 2);
   ExpectSameContainment(original, *optimized);
}

TEST(TGeoCompositeShapeOptimize, RejectsIntersection)
{
   TGeoBBox first(2., 2., 2.);
   TGeoBBox second(2., 2., 2.);
   TGeoBBox third(2., 2., 2.);
   TGeoCompositeShape nested("intersection_union", new TGeoUnion(&first, &second));
   TGeoCompositeShape original("intersection", new TGeoIntersection(&nested, &third));

   EXPECT_EQ(original.Optimize(), &original);
}

TEST(TGeoCompositeShapeOptimize, RejectsInfiniteHalfSpaceComponent)
{
   TGeoBBox outer(3., 3., 3.);
   TGeoBBox hole(1., 1., 4.);
   double planePoint[3] = {0., 0., 0.};
   double planeNormal[3] = {1., 0., 0.};
   TGeoHalfSpace halfSpace("clipping_plane", planePoint, planeNormal);
   TGeoCompositeShape clipped("clipped", new TGeoSubtraction(&outer, &halfSpace));
   TGeoCompositeShape original("clipped_with_hole", new TGeoSubtraction(&clipped, &hole));

   EXPECT_FALSE(original.CanOptimize());
   EXPECT_EQ(original.Optimize(), &original);
}

TEST(TGeoCompositeShapeOptimize, RejectsSubtractionOnRightHandSide)
{
   TGeoBBox first(3., 3., 3.);
   TGeoBBox second(2., 2., 2.);
   TGeoBBox third(1., 1., 1.);
   TGeoCompositeShape nested("nested_subtraction", new TGeoSubtraction(&second, &third));
   TGeoCompositeShape original("right_subtraction", new TGeoSubtraction(&first, &nested));

   EXPECT_EQ(original.Optimize(), &original);
}
