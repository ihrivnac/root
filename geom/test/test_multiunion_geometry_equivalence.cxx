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
#include <TGeoNavigator.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TObjArray.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using Vec3 = std::array<Double_t, 3>;
using PathId = std::uint32_t;

constexpr PathId kInvalidPath = std::numeric_limits<PathId>::max();
constexpr Double_t kTrackLimit = 1.e5;
constexpr Double_t kDistanceAbsTolerance = 5.e-6;
constexpr Double_t kDistanceRelTolerance = 1.e-9;

struct TestConfig {
   std::uint64_t fSeed{0x5eed1234ULL};
   std::size_t fDepth{4};
   std::size_t fPointCount{2000};
   std::size_t fRayCount{500};
   std::size_t fMaxBoundaries{64};
};

struct PointQuery {
   Vec3 fPoint;
};

struct RayQuery {
   Vec3 fPoint;
   Vec3 fDirection;
};

struct QueryCorpus {
   std::vector<PointQuery> fPoints;
   std::vector<RayQuery> fRays;
};

struct PointResponse {
   PathId fPath{kInvalidPath};
   Double_t fSafety{0.};
};

enum class RayStatus : unsigned char {
   kExited,
   kNoBoundary,
   kIterationLimit,
   kRepeatedSmallSteps
};

struct BoundaryResponse {
   PathId fPathBefore{kInvalidPath};
   PathId fPathAfter{kInvalidPath};
   Double_t fDistance{0.};
   Vec3 fPoint{};
};

struct RayResponse {
   RayStatus fStatus{RayStatus::kNoBoundary};
   std::vector<BoundaryResponse> fBoundaries;
};

struct GeometrySnapshot {
   std::vector<PointResponse> fPoints;
   std::vector<RayResponse> fRays;
};

struct TimingResult {
   Double_t fPointMilliseconds{0.};
   Double_t fRayMilliseconds{0.};
   Double_t fSafetyChecksum{0.};
   Double_t fStepChecksum{0.};
   std::size_t fCrossings{0};
};

struct NormalResponse {
   Vec3 fNormal{};
};

struct ShapeQueryCorpus {
   std::vector<PointQuery> fPoints;
   std::vector<RayQuery> fRays;
};

struct ShapeResponse {
   Bool_t fInside{kFALSE};
   Double_t fSafety{0.};
   Double_t fDistance{0.};
};

struct ShapeTimingResult {
   Double_t fPointMilliseconds{0.};
   Double_t fRayMilliseconds{0.};
   Double_t fChecksum{0.};
};

struct ShapeCandidate {
   TGeoCompositeShape *fOriginal{nullptr};
   TGeoShape *fOptimized{nullptr};
   std::size_t fLeaves{0};
   std::string fName;
};

struct ShapeTimingBin {
   std::size_t fMinimumLeaves{0};
   std::size_t fMaximumLeaves{0};
   std::size_t fShapeCount{0};
   Double_t fOriginalPointMilliseconds{0.};
   Double_t fOptimizedPointMilliseconds{0.};
   Double_t fOriginalRayMilliseconds{0.};
   Double_t fOptimizedRayMilliseconds{0.};
};

class PathTable {
private:
   std::vector<std::string> fPaths;
   std::unordered_map<std::string, PathId> fIds;

public:
   PathId Intern(const std::string &path)
   {
      auto found = fIds.find(path);
      if (found != fIds.end())
         return found->second;
      const auto id = static_cast<PathId>(fPaths.size());
      fPaths.push_back(path);
      fIds.emplace(fPaths.back(), id);
      return id;
   }

   PathId Find(const std::string &path) const
   {
      auto found = fIds.find(path);
      return found == fIds.end() ? kInvalidPath : found->second;
   }

   const std::string &Get(PathId id) const
   {
      static const std::string unknown{"<unknown>"};
      return id < fPaths.size() ? fPaths[id] : unknown;
   }
};

struct GeometryFixture {
   TGeoVolume *fTop{nullptr};
   TGeoVolume *fSharedFirst{nullptr};
   TGeoVolume *fSharedSecond{nullptr};
   TGeoVolume *fShapeMinusUnion{nullptr};
   TGeoVolume *fUnionMinusUnion{nullptr};
   TGeoVolume *fIneligible{nullptr};
   TGeoVolume *fUnused{nullptr};
   TGeoCompositeShape *fSharedOriginal{nullptr};
   TGeoCompositeShape *fShapeMinusUnionOriginal{nullptr};
   TGeoCompositeShape *fUnionMinusUnionOriginal{nullptr};
   TGeoCompositeShape *fIneligibleOriginal{nullptr};
   TGeoCompositeShape *fUnusedOriginal{nullptr};
};

std::size_t ReadSize(const char *name, std::size_t fallback)
{
   const char *text = std::getenv(name);
   if (!text || !text[0])
      return fallback;
   char *end = nullptr;
   const auto value = std::strtoull(text, &end, 10);
   return end != text && *end == '\0' && value > 0 ? static_cast<std::size_t>(value) : fallback;
}

std::uint64_t ReadSeed(const char *name, std::uint64_t fallback)
{
   const char *text = std::getenv(name);
   if (!text || !text[0])
      return fallback;
   char *end = nullptr;
   const auto value = std::strtoull(text, &end, 0);
   return end != text && *end == '\0' ? static_cast<std::uint64_t>(value) : fallback;
}

TestConfig GetConfig()
{
   TestConfig config;
   const char *stress = std::getenv("TGEOMU_STRESS");
   if (stress && stress[0] && std::string(stress) != "0") {
      config.fPointCount = 100000;
      config.fRayCount = 10000;
      config.fMaxBoundaries = 128;
   }
   config.fSeed = ReadSeed("TGEOMU_SEED", config.fSeed);
   config.fDepth = ReadSize("TGEOMU_DEPTH", config.fDepth);
   config.fPointCount = ReadSize("TGEOMU_POINTS", config.fPointCount);
   config.fRayCount = ReadSize("TGEOMU_RAYS", config.fRayCount);
   config.fMaxBoundaries = ReadSize("TGEOMU_MAX_BOUNDARIES", config.fMaxBoundaries);
   return config;
}

Bool_t TimingEnabled()
{
   const char *timing = std::getenv("TGEOMU_TIMING");
   return timing && timing[0] && std::string(timing) != "0";
}

const char *GeometryFile()
{
   const char *filename = std::getenv("TGEOMU_GEOMETRY_FILE");
   return filename && filename[0] ? filename : nullptr;
}

TGeoTranslation *Translation(Double_t x, Double_t y = 0., Double_t z = 0.)
{
   return new TGeoTranslation(x, y, z);
}

Vec3 ComponentCenter(std::size_t index)
{
   constexpr std::size_t columns = 10;
   return {-5.4 + 1.2 * static_cast<Double_t>(index % columns),
           -2.4 + 1.2 * static_cast<Double_t>(index / columns), 0.};
}

TGeoCompositeShape *MakeBoxUnion(const char *prefix, std::size_t depth, const Vec3 &halfLengths)
{
   const auto makeBox = [&](std::size_t index) {
      const std::string name = std::string(prefix) + "_box_" + std::to_string(index);
      return new TGeoBBox(name.c_str(), halfLengths[0], halfLengths[1], halfLengths[2]);
   };
   const auto matrix = [](std::size_t index) {
      const Vec3 center = ComponentCenter(index);
      return Translation(center[0], center[1], center[2]);
   };

   auto *first = makeBox(0);
   auto *second = makeBox(1);
   std::string name = std::string(prefix) + "_union_2";
   TGeoCompositeShape *result =
      new TGeoCompositeShape(name.c_str(), new TGeoUnion(first, second, matrix(0), matrix(1)));
   for (std::size_t index = 2; index < depth; ++index) {
      name = std::string(prefix) + "_union_" + std::to_string(index + 1);
      result = new TGeoCompositeShape(name.c_str(), new TGeoUnion(result, makeBox(index), nullptr, matrix(index)));
   }
   return result;
}

TGeoCompositeShape *MakeShapeMinusUnion(std::size_t depth)
{
   auto *outer = new TGeoBBox("equiv_sub_outer", 7., 5., 5.);
   auto *holes = MakeBoxUnion("equiv_sub_hole", depth, {.22, .22, 6.});
   return new TGeoCompositeShape("equiv_shape_minus_union", new TGeoSubtraction(outer, holes));
}

TGeoCompositeShape *MakeUnionMinusUnion(std::size_t depth)
{
   auto *positive = MakeBoxUnion("equiv_positive", depth, {.5, .5, .75});
   auto *negative = MakeBoxUnion("equiv_negative", depth, {.16, .16, 1.});
   return new TGeoCompositeShape("equiv_union_minus_union", new TGeoSubtraction(positive, negative));
}

TGeoVolume *MakeContainer(const char *name, TGeoMedium *medium, TGeoVolume *content)
{
   const std::string shapeName = std::string(name) + "_shape";
   auto *container = new TGeoVolume(name, new TGeoBBox(shapeName.c_str(), 12., 9., 9.), medium);
   container->AddNode(content, 1);
   return container;
}

GeometryFixture BuildGeometry(TGeoManager &manager, const TestConfig &config)
{
   GeometryFixture fixture;
   auto *material = new TGeoMaterial("equiv_material", 1., 1., 1.);
   auto *medium = new TGeoMedium("equiv_medium", 1, material);
   fixture.fTop = new TGeoVolume("equiv_top", new TGeoBBox("equiv_top_shape", 100., 40., 30.), medium);
   manager.SetTopVolume(fixture.fTop);

   fixture.fSharedOriginal = MakeBoxUnion("equiv_union", config.fDepth, {.5, .5, .5});
   fixture.fSharedFirst = new TGeoVolume("equiv_shared_first", fixture.fSharedOriginal, medium);
   fixture.fSharedSecond = new TGeoVolume("equiv_shared_second", fixture.fSharedOriginal, medium);
   fixture.fShapeMinusUnionOriginal = MakeShapeMinusUnion(config.fDepth);
   fixture.fShapeMinusUnion =
      new TGeoVolume("equiv_shape_minus_union_volume", fixture.fShapeMinusUnionOriginal, medium);
   fixture.fUnionMinusUnionOriginal = MakeUnionMinusUnion(config.fDepth);
   fixture.fUnionMinusUnion =
      new TGeoVolume("equiv_union_minus_union_volume", fixture.fUnionMinusUnionOriginal, medium);

   auto *ineligibleLeft = new TGeoBBox("equiv_ineligible_left", 3., 2., 2.);
   auto *ineligibleRight = new TGeoBBox("equiv_ineligible_right", 2., 3., 2.);
   fixture.fIneligibleOriginal = new TGeoCompositeShape(
      "equiv_ineligible_intersection", new TGeoIntersection(ineligibleLeft, ineligibleRight));
   fixture.fIneligible = new TGeoVolume("equiv_ineligible_volume", fixture.fIneligibleOriginal, medium);

   auto *unusedFirst = new TGeoBBox("equiv_unused_first", 1., 1., 1.);
   auto *unusedSecond = new TGeoBBox("equiv_unused_second", 1., 1., 1.);
   auto *unusedThird = new TGeoBBox("equiv_unused_third", 1., 1., 1.);
   auto *unused12 = new TGeoCompositeShape("equiv_unused12", new TGeoUnion(unusedFirst, unusedSecond));
   fixture.fUnusedOriginal = new TGeoCompositeShape("equiv_unused", new TGeoUnion(unused12, unusedThird));
   fixture.fUnused = new TGeoVolume("equiv_unused_volume", fixture.fUnusedOriginal, medium);

   fixture.fTop->AddNode(MakeContainer("equiv_container_shared_first", medium, fixture.fSharedFirst), 1,
                         Translation(-60.));
   fixture.fTop->AddNode(MakeContainer("equiv_container_shared_second", medium, fixture.fSharedSecond), 1,
                         Translation(-20.));
   fixture.fTop->AddNode(MakeContainer("equiv_container_subtraction", medium, fixture.fShapeMinusUnion), 1,
                         Translation(20.));
   fixture.fTop->AddNode(MakeContainer("equiv_container_difference", medium, fixture.fUnionMinusUnion), 1,
                         Translation(60.));
   fixture.fTop->AddNode(MakeContainer("equiv_container_ineligible", medium, fixture.fIneligible), 1,
                         Translation(0., 22.));
   manager.CloseGeometry();
   return fixture;
}

Double_t Uniform(std::mt19937_64 &engine, Double_t low, Double_t high)
{
   return std::uniform_real_distribution<Double_t>(low, high)(engine);
}

Vec3 RandomDirection(std::mt19937_64 &engine)
{
   const Double_t z = Uniform(engine, -1., 1.);
   const Double_t phi = Uniform(engine, 0., 2. * std::acos(-1.));
   const Double_t radial = std::sqrt(std::max(0., 1. - z * z));
   return {radial * std::cos(phi), radial * std::sin(phi), z};
}

QueryCorpus MakeQueries(const TestConfig &config)
{
   QueryCorpus corpus;
   corpus.fPoints.reserve(config.fPointCount);
   corpus.fRays.reserve(config.fRayCount);
   std::mt19937_64 engine(config.fSeed);
   constexpr Double_t centers[] = {-60., -20., 20., 60.};

   // Put deterministic probes immediately on both sides of component faces.
   for (Double_t center : centers) {
      for (std::size_t component = 0; component < config.fDepth; ++component) {
         const Vec3 local = ComponentCenter(component);
         for (Double_t halfLength : {.5, .22, .16}) {
            for (Double_t side : {-1., 1.}) {
               for (Double_t offset : {-1.e-7, 1.e-7}) {
                  if (corpus.fPoints.size() < config.fPointCount)
                     corpus.fPoints.push_back(
                        {{center + local[0] + side * halfLength + offset, local[1], local[2]}});
               }
            }
         }
      }
   }

   while (corpus.fPoints.size() < config.fPointCount) {
      if ((corpus.fPoints.size() & 1U) == 0U) {
         const Double_t center = centers[corpus.fPoints.size() % 4];
         corpus.fPoints.push_back(
            {{Uniform(engine, center - 8., center + 8.), Uniform(engine, -6., 6.), Uniform(engine, -6., 6.)}});
      } else {
         corpus.fPoints.push_back(
            {{Uniform(engine, -110., 110.), Uniform(engine, -45., 45.), Uniform(engine, -35., 35.)}});
      }
   }

   // Focused tracks cross complete component rows in every placed optimizable solid.
   const std::size_t rows = (config.fDepth + 9) / 10;
   for (Double_t center : centers) {
      for (std::size_t row = 0; row < rows && corpus.fRays.size() < config.fRayCount; ++row) {
         const Double_t y = ComponentCenter(10 * row).at(1);
         corpus.fRays.push_back({{center - 11., y, 0.}, {1., 0., 0.}});
      }
   }

   while (corpus.fRays.size() < config.fRayCount) {
      const std::size_t index = corpus.fRays.size();
      if (index % 4 == 0) {
         corpus.fRays.push_back({{-110., Uniform(engine, -8., 8.), Uniform(engine, -8., 8.)}, {1., 0., 0.}});
      } else if (index % 4 == 1) {
         const Double_t center = centers[(index / 4) % 4];
         corpus.fRays.push_back(
            {{Uniform(engine, center - 7., center + 7.), Uniform(engine, -5., 5.), Uniform(engine, -5., 5.)},
             RandomDirection(engine)});
      } else {
         corpus.fRays.push_back(
            {{Uniform(engine, -95., 95.), Uniform(engine, -35., 35.), Uniform(engine, -25., 25.)},
             RandomDirection(engine)});
      }
   }
   return corpus;
}

std::string CurrentPath(TGeoManager &manager)
{
   if (manager.IsOutside() || !manager.GetCurrentNode())
      return "<outside>";
   return manager.GetPath();
}

PathId ResolvePath(PathTable &paths, const std::string &path, Bool_t addPaths)
{
   return addPaths ? paths.Intern(path) : paths.Find(path);
}

Bool_t IsNear(Double_t first, Double_t second)
{
   return std::abs(first - second) <=
          kDistanceAbsTolerance + kDistanceRelTolerance * std::max(std::abs(first), std::abs(second));
}

PointResponse CapturePoint(TGeoManager &manager, const PointQuery &query, PathTable &paths, Bool_t addPaths)
{
   manager.FindNode(query.fPoint[0], query.fPoint[1], query.fPoint[2]);
   PointResponse response;
   response.fPath = ResolvePath(paths, CurrentPath(manager), addPaths);
   response.fSafety = manager.Safety();
   return response;
}

void ValidateSafety(TGeoManager &manager, const PointQuery &query, const PointResponse &response, const PathTable &paths,
                    std::size_t queryId, std::uint64_t seed)
{
   SCOPED_TRACE(::testing::Message() << "seed=" << seed << " point query=" << queryId << " safety="
                                     << response.fSafety);
   EXPECT_TRUE(std::isfinite(response.fSafety));
   EXPECT_GE(response.fSafety, 0.);
   if (response.fPath == kInvalidPath || response.fSafety <= 1.e-8 || response.fSafety >= TGeoShape::Big())
      return;

   const Double_t scale = .5 * response.fSafety;
   const Vec3 direction = {0.2672612419124244, -0.5345224838248488, 0.8017837257372732};
   for (Double_t sign : {-1., 1.}) {
      Vec3 probe = {query.fPoint[0] + sign * scale * direction[0],
                    query.fPoint[1] + sign * scale * direction[1],
                    query.fPoint[2] + sign * scale * direction[2]};
      manager.FindNode(probe[0], probe[1], probe[2]);
      EXPECT_EQ(CurrentPath(manager), paths.Get(response.fPath));
   }
}

RayResponse CaptureRay(TGeoManager &manager, const RayQuery &query, std::size_t maxBoundaries, PathTable &paths,
                       Bool_t addPaths)
{
   RayResponse response;
   response.fBoundaries.reserve(16);
   manager.InitTrack(query.fPoint.data(), query.fDirection.data());
   Int_t smallSteps = 0;

   for (std::size_t crossing = 0; crossing < maxBoundaries; ++crossing) {
      BoundaryResponse boundary;
      boundary.fPathBefore = ResolvePath(paths, CurrentPath(manager), addPaths);
      const Vec3 start = {manager.GetCurrentPoint()[0], manager.GetCurrentPoint()[1], manager.GetCurrentPoint()[2]};

      TGeoNode *nextBoundary = manager.FindNextBoundary(kTrackLimit);
      const Double_t distance = manager.GetStep();
      if (!nextBoundary || distance >= .999 * kTrackLimit || distance >= TGeoShape::Big()) {
         response.fStatus = RayStatus::kNoBoundary;
         return response;
      }

      boundary.fDistance = distance;
      boundary.fPoint = {start[0] + distance * query.fDirection[0], start[1] + distance * query.fDirection[1],
                         start[2] + distance * query.fDirection[2]};

      manager.Step(kTRUE, kTRUE);
      boundary.fPathAfter = ResolvePath(paths, CurrentPath(manager), addPaths);
      response.fBoundaries.push_back(boundary);

      if (distance < 1.e-8)
         ++smallSteps;
      else
         smallSteps = 0;
      if (smallSteps > 3) {
         response.fStatus = RayStatus::kRepeatedSmallSteps;
         return response;
      }
      if (manager.IsOutside()) {
         response.fStatus = RayStatus::kExited;
         return response;
      }
   }
   response.fStatus = RayStatus::kIterationLimit;
   return response;
}

GeometrySnapshot CaptureBaseline(TGeoManager &manager, const QueryCorpus &corpus, const TestConfig &config,
                                 PathTable &paths)
{
   GeometrySnapshot snapshot;
   snapshot.fPoints.reserve(corpus.fPoints.size());
   snapshot.fRays.reserve(corpus.fRays.size());
   paths.Intern("<outside>");

   for (std::size_t index = 0; index < corpus.fPoints.size(); ++index) {
      auto response = CapturePoint(manager, corpus.fPoints[index], paths, kTRUE);
      ValidateSafety(manager, corpus.fPoints[index], response, paths, index, config.fSeed);
      snapshot.fPoints.push_back(response);
   }
   for (const auto &query : corpus.fRays)
      snapshot.fRays.push_back(CaptureRay(manager, query, config.fMaxBoundaries, paths, kTRUE));
   return snapshot;
}

void ComparePointResponses(TGeoManager &manager, const QueryCorpus &corpus, const GeometrySnapshot &baseline,
                           const TestConfig &config, PathTable &paths)
{
   ASSERT_EQ(corpus.fPoints.size(), baseline.fPoints.size());
   for (std::size_t index = 0; index < corpus.fPoints.size(); ++index) {
      SCOPED_TRACE(::testing::Message() << "seed=" << config.fSeed << " point query=" << index << " point=("
                                        << corpus.fPoints[index].fPoint[0] << ", " << corpus.fPoints[index].fPoint[1]
                                        << ", " << corpus.fPoints[index].fPoint[2] << ")");
      auto optimized = CapturePoint(manager, corpus.fPoints[index], paths, kFALSE);
      EXPECT_NE(optimized.fPath, kInvalidPath);
      EXPECT_EQ(optimized.fPath, baseline.fPoints[index].fPath)
         << "original=" << paths.Get(baseline.fPoints[index].fPath)
         << " optimized=" << CurrentPath(manager);
      ValidateSafety(manager, corpus.fPoints[index], optimized, paths, index, config.fSeed);
   }
}

void CompareRayResponses(TGeoManager &manager, const QueryCorpus &corpus, const GeometrySnapshot &baseline,
                         const TestConfig &config, PathTable &paths)
{
   ASSERT_EQ(corpus.fRays.size(), baseline.fRays.size());
   for (std::size_t ray = 0; ray < corpus.fRays.size(); ++ray) {
      SCOPED_TRACE(::testing::Message() << "seed=" << config.fSeed << " ray=" << ray << " point=("
                                        << corpus.fRays[ray].fPoint[0] << ", " << corpus.fRays[ray].fPoint[1] << ", "
                                        << corpus.fRays[ray].fPoint[2] << ") direction=("
                                        << corpus.fRays[ray].fDirection[0] << ", "
                                        << corpus.fRays[ray].fDirection[1] << ", "
                                        << corpus.fRays[ray].fDirection[2] << ")");
      const auto optimized = CaptureRay(manager, corpus.fRays[ray], config.fMaxBoundaries, paths, kFALSE);
      const auto &original = baseline.fRays[ray];
      EXPECT_EQ(optimized.fStatus, original.fStatus);
      ASSERT_EQ(optimized.fBoundaries.size(), original.fBoundaries.size());

      for (std::size_t crossing = 0; crossing < original.fBoundaries.size(); ++crossing) {
         SCOPED_TRACE(::testing::Message() << "boundary=" << crossing);
         const auto &first = original.fBoundaries[crossing];
         const auto &second = optimized.fBoundaries[crossing];
         EXPECT_NE(second.fPathBefore, kInvalidPath);
         EXPECT_NE(second.fPathAfter, kInvalidPath);
         EXPECT_EQ(second.fPathBefore, first.fPathBefore)
            << "original=" << paths.Get(first.fPathBefore) << " optimized=" << paths.Get(second.fPathBefore);
         EXPECT_EQ(second.fPathAfter, first.fPathAfter)
            << "original=" << paths.Get(first.fPathAfter) << " optimized=" << paths.Get(second.fPathAfter);
         EXPECT_TRUE(IsNear(second.fDistance, first.fDistance))
            << "original distance=" << first.fDistance << " optimized distance=" << second.fDistance;
         for (Int_t axis = 0; axis < 3; ++axis)
            EXPECT_TRUE(IsNear(second.fPoint[axis], first.fPoint[axis]));

      }
   }
}

TimingResult MeasureNavigation(TGeoManager &manager, const QueryCorpus &corpus, const TestConfig &config)
{
   using Clock = std::chrono::steady_clock;
   TimingResult result;

   auto start = Clock::now();
   for (const auto &query : corpus.fPoints) {
      manager.FindNode(query.fPoint[0], query.fPoint[1], query.fPoint[2]);
      result.fSafetyChecksum += manager.Safety();
   }
   auto stop = Clock::now();
   result.fPointMilliseconds = std::chrono::duration<Double_t, std::milli>(stop - start).count();

   start = Clock::now();
   for (const auto &query : corpus.fRays) {
      manager.InitTrack(query.fPoint.data(), query.fDirection.data());
      Int_t smallSteps = 0;
      for (std::size_t crossing = 0; crossing < config.fMaxBoundaries; ++crossing) {
         TGeoNode *nextBoundary = manager.FindNextBoundary(kTrackLimit);
         const Double_t distance = manager.GetStep();
         if (!nextBoundary || distance >= .999 * kTrackLimit || distance >= TGeoShape::Big())
            break;
         result.fStepChecksum += distance;
         ++result.fCrossings;
         manager.Step(kTRUE, kTRUE);
         if (distance < 1.e-8)
            ++smallSteps;
         else
            smallSteps = 0;
         if (smallSteps > 3 || manager.IsOutside())
            break;
      }
   }
   stop = Clock::now();
   result.fRayMilliseconds = std::chrono::duration<Double_t, std::milli>(stop - start).count();
   return result;
}

void PrintTimingSummary(const TestConfig &config, const TimingResult &original, const TimingResult &optimized)
{
   const auto speedup = [](Double_t before, Double_t after) { return after > 0. ? before / after : 0.; };
   std::cout << std::fixed << std::setprecision(3)
             << "\nTGeoMultiUnion navigation timing (depth=" << config.fDepth << ", seed=" << config.fSeed
             << ", points=" << config.fPointCount << ", rays=" << config.fRayCount << ")\n"
             << "  Original geometry:  point location+safety " << original.fPointMilliseconds
             << " ms, ray navigation " << original.fRayMilliseconds << " ms\n"
             << "  Optimized geometry: point location+safety " << optimized.fPointMilliseconds
             << " ms, ray navigation " << optimized.fRayMilliseconds << " ms\n"
             << "  Speedup:             point location+safety "
             << speedup(original.fPointMilliseconds, optimized.fPointMilliseconds) << "x, ray navigation "
             << speedup(original.fRayMilliseconds, optimized.fRayMilliseconds) << "x\n"
             << "  Navigation crossings: original " << original.fCrossings << ", optimized "
             << optimized.fCrossings << "\n";
}

std::size_t CountPrimitiveLeaves(const TGeoShape &shape)
{
   auto *composite = dynamic_cast<const TGeoCompositeShape *>(&shape);
   if (!composite || !composite->GetBoolNode())
      return 1;
   return CountPrimitiveLeaves(*composite->GetBoolNode()->GetLeftShape()) +
          CountPrimitiveLeaves(*composite->GetBoolNode()->GetRightShape());
}

std::vector<ShapeCandidate> FindShapeCandidates(TGeoManager &manager, std::size_t minimumDepth)
{
   std::vector<ShapeCandidate> candidates;
   std::unordered_set<TGeoCompositeShape *> seen;
   TObjArray *shapes = manager.GetListOfShapes();
   if (!shapes)
      return candidates;

   const Int_t entries = shapes->GetEntriesFast();
   for (Int_t index = 0; index < entries; ++index) {
      auto *composite = dynamic_cast<TGeoCompositeShape *>(shapes->At(index));
      if (!composite || !seen.insert(composite).second || !composite->CanOptimize())
         continue;
      const std::size_t leaves = CountPrimitiveLeaves(*composite);
      if (leaves < minimumDepth)
         continue;
      const char *name = composite->GetName();
      candidates.push_back({composite, nullptr, leaves, name && name[0] ? name : "<unnamed>"});
   }
   return candidates;
}

ShapeQueryCorpus MakeShapeQueries(const TGeoCompositeShape &shape, std::size_t pointCount, std::size_t rayCount,
                                  std::uint64_t seed)
{
   ShapeQueryCorpus corpus;
   corpus.fPoints.reserve(pointCount);
   corpus.fRays.reserve(rayCount);
   std::mt19937_64 engine(seed);
   const Double_t *origin = shape.GetOrigin();
   const Vec3 halfLengths = {std::max(shape.GetDX(), 1.e-9), std::max(shape.GetDY(), 1.e-9),
                             std::max(shape.GetDZ(), 1.e-9)};

   for (Int_t axis = 0; axis < 3 && corpus.fPoints.size() < pointCount; ++axis) {
      for (Double_t side : {-1., 1.}) {
         for (Double_t offset : {-1.e-7, 1.e-7}) {
            Vec3 point{origin[0], origin[1], origin[2]};
            point[axis] += side * halfLengths[axis] + offset * std::max(1., halfLengths[axis]);
            corpus.fPoints.push_back({point});
         }
      }
   }
   while (corpus.fPoints.size() < pointCount) {
      const Double_t scale = (corpus.fPoints.size() & 1U) ? 1.15 : .999999;
      corpus.fPoints.push_back({{origin[0] + Uniform(engine, -scale, scale) * halfLengths[0],
                                 origin[1] + Uniform(engine, -scale, scale) * halfLengths[1],
                                 origin[2] + Uniform(engine, -scale, scale) * halfLengths[2]}});
   }
   while (corpus.fRays.size() < rayCount) {
      const Double_t scale = (corpus.fRays.size() & 1U) ? 1.2 : .95;
      corpus.fRays.push_back({{origin[0] + Uniform(engine, -scale, scale) * halfLengths[0],
                               origin[1] + Uniform(engine, -scale, scale) * halfLengths[1],
                               origin[2] + Uniform(engine, -scale, scale) * halfLengths[2]},
                              RandomDirection(engine)});
   }
   return corpus;
}

ShapeResponse CaptureShapeResponse(TGeoShape &shape, const Vec3 &point, const Vec3 &direction)
{
   ShapeResponse response;
   response.fInside = shape.Contains(point.data());
   response.fSafety = shape.Safety(point.data(), response.fInside);
   response.fDistance = response.fInside ? shape.DistFromInside(point.data(), direction.data(), 3)
                                         : shape.DistFromOutside(point.data(), direction.data(), 3);
   return response;
}

void ValidateShapeSafety(TGeoShape &shape, const Vec3 &point, const ShapeResponse &response)
{
   EXPECT_TRUE(std::isfinite(response.fSafety));
   // Some existing ROOT primitives return a negative safety estimate for points
   // close to, or numerically classified on the wrong side of, their surface.
   // This is independent of the composite rewrite, so only validate the
   // containment guarantee when a usable positive safety was returned.
   if (response.fSafety <= 1.e-9 || response.fSafety >= TGeoShape::Big())
      return;
   constexpr Vec3 direction{0.2672612419124244, -0.5345224838248488, 0.8017837257372732};
   const Double_t distance = .5 * response.fSafety;
   for (Double_t sign : {-1., 1.}) {
      const Vec3 probe{point[0] + sign * distance * direction[0], point[1] + sign * distance * direction[1],
                       point[2] + sign * distance * direction[2]};
      EXPECT_EQ(shape.Contains(probe.data()), response.fInside);
   }
}

ShapeTimingResult MeasureShape(TGeoShape &shape, const ShapeQueryCorpus &corpus)
{
   using Clock = std::chrono::steady_clock;
   ShapeTimingResult result;
   auto start = Clock::now();
   for (const auto &query : corpus.fPoints) {
      const Bool_t inside = shape.Contains(query.fPoint.data());
      result.fChecksum += shape.Safety(query.fPoint.data(), inside);
   }
   auto stop = Clock::now();
   result.fPointMilliseconds = std::chrono::duration<Double_t, std::milli>(stop - start).count();

   start = Clock::now();
   for (const auto &query : corpus.fRays) {
      const Bool_t inside = shape.Contains(query.fPoint.data());
      result.fChecksum += inside ? shape.DistFromInside(query.fPoint.data(), query.fDirection.data(), 3)
                                 : shape.DistFromOutside(query.fPoint.data(), query.fDirection.data(), 3);
   }
   stop = Clock::now();
   result.fRayMilliseconds = std::chrono::duration<Double_t, std::milli>(stop - start).count();
   return result;
}

void CompareShape(TGeoShape &original, TGeoShape &optimized, const ShapeQueryCorpus &corpus,
                  const std::string &name, std::uint64_t seed)
{
   for (std::size_t index = 0; index < corpus.fPoints.size(); ++index) {
      SCOPED_TRACE(::testing::Message() << "shape=" << name << " seed=" << seed << " point=" << index);
      const Vec3 &point = corpus.fPoints[index].fPoint;
      const Vec3 direction = corpus.fRays[index % corpus.fRays.size()].fDirection;
      const auto before = CaptureShapeResponse(original, point, direction);
      const auto after = CaptureShapeResponse(optimized, point, direction);
      EXPECT_EQ(after.fInside, before.fInside);
      EXPECT_TRUE(IsNear(after.fDistance, before.fDistance))
         << "original distance=" << before.fDistance << " optimized distance=" << after.fDistance;
      ValidateShapeSafety(original, point, before);
      ValidateShapeSafety(optimized, point, after);
   }
   for (std::size_t index = 0; index < corpus.fRays.size(); ++index) {
      SCOPED_TRACE(::testing::Message() << "shape=" << name << " seed=" << seed << " ray=" << index);
      const auto &query = corpus.fRays[index];
      const auto before = CaptureShapeResponse(original, query.fPoint, query.fDirection);
      const auto after = CaptureShapeResponse(optimized, query.fPoint, query.fDirection);
      EXPECT_EQ(after.fInside, before.fInside);
      EXPECT_TRUE(IsNear(after.fDistance, before.fDistance))
         << "original distance=" << before.fDistance << " optimized distance=" << after.fDistance;
   }
}

std::array<ShapeTimingBin, 5> MakeShapeTimingBins()
{
   return {{{3, 5}, {6, 9}, {10, 19}, {20, 49}, {50, std::numeric_limits<std::size_t>::max()}}};
}

ShapeTimingBin &FindShapeTimingBin(std::array<ShapeTimingBin, 5> &bins, std::size_t leaves)
{
   for (auto &bin : bins) {
      if (leaves >= bin.fMinimumLeaves && leaves <= bin.fMaximumLeaves)
         return bin;
   }
   return bins.back();
}

void PrintShapeTimingBins(const std::array<ShapeTimingBin, 5> &bins, std::size_t pointsPerShape,
                          std::size_t raysPerShape)
{
   const auto speedup = [](Double_t original, Double_t optimized) {
      return optimized > 0. ? original / optimized : 0.;
   };
   std::cout << "\nTGeoMultiUnion imported geometry timing by primitive-leaf count\n"
             << "  Queries per shape: " << pointsPerShape << " points, " << raysPerShape << " rays\n"
             << "  Leaves    Shapes       Point time original -> optimized (speedup)"
                "       Ray time original -> optimized (speedup)\n";
   for (const auto &bin : bins) {
      if (!bin.fShapeCount)
         continue;
      const std::string label =
         bin.fMaximumLeaves == std::numeric_limits<std::size_t>::max()
            ? std::to_string(bin.fMinimumLeaves) + "+"
            : std::to_string(bin.fMinimumLeaves) + "-" + std::to_string(bin.fMaximumLeaves);
      std::cout << "  " << std::setw(7) << std::left << label << std::right << std::setw(8) << bin.fShapeCount
                << std::fixed << std::setprecision(3) << std::setw(16) << bin.fOriginalPointMilliseconds << " -> "
                << std::setw(10) << bin.fOptimizedPointMilliseconds << " ms ("
                << std::setw(7) << speedup(bin.fOriginalPointMilliseconds, bin.fOptimizedPointMilliseconds)
                << "x)" << std::setw(16) << bin.fOriginalRayMilliseconds << " -> " << std::setw(10)
                << bin.fOptimizedRayMilliseconds << " ms ("
                << std::setw(7) << speedup(bin.fOriginalRayMilliseconds, bin.fOptimizedRayMilliseconds) << "x)\n";
   }
}

void RunImportedGeometryTest(const TestConfig &config, const char *filename)
{
   const char *keyName = std::getenv("TGEOMU_GEOMETRY_NAME");
   TGeoManager *manager = TGeoManager::Import(filename, keyName ? keyName : "");
   ASSERT_NE(manager, nullptr) << "Could not import geometry file " << filename;

   auto candidates = FindShapeCandidates(*manager, config.fDepth);
   std::cout << "\nTGeoMultiUnion imported geometry scan\n"
             << "  File: " << filename << "\n"
             << "  Minimum component depth: " << config.fDepth << "\n"
             << "  Optimizable composite shapes selected: " << candidates.size() << "\n";
   for (std::size_t index = 0; index < candidates.size(); ++index)
      std::cout << "    [" << index << "] " << candidates[index].fName << " (leaves=" << candidates[index].fLeaves
                << ")\n";
   if (candidates.empty()) {
      GTEST_SKIP() << "No optimizable composite shapes satisfy TGEOMU_DEPTH=" << config.fDepth;
   }

   const std::size_t pointsPerShape = std::max<std::size_t>(256, config.fPointCount / candidates.size());
   const std::size_t raysPerShape = std::max<std::size_t>(64, config.fRayCount / candidates.size());
   const Bool_t timingEnabled = TimingEnabled();
   Double_t originalPointTime = 0.;
   Double_t optimizedPointTime = 0.;
   Double_t originalRayTime = 0.;
   Double_t optimizedRayTime = 0.;
   auto timingBins = MakeShapeTimingBins();

   for (std::size_t index = 0; index < candidates.size(); ++index) {
      auto &candidate = candidates[index];
      const std::uint64_t seed = config.fSeed ^ (0x9e3779b97f4a7c15ULL * (index + 1));
      const ShapeQueryCorpus corpus = MakeShapeQueries(*candidate.fOriginal, pointsPerShape, raysPerShape, seed);
      candidate.fOptimized = candidate.fOriginal->Optimize();
      ASSERT_NE(candidate.fOptimized, candidate.fOriginal);
      const std::string diagnosticName = "[" + std::to_string(index) + "] " + candidate.fName +
                                         " (leaves=" + std::to_string(candidate.fLeaves) + ")";
      CompareShape(*candidate.fOriginal, *candidate.fOptimized, corpus, diagnosticName, seed);
      if (timingEnabled) {
         const auto before = MeasureShape(*candidate.fOriginal, corpus);
         const auto after = MeasureShape(*candidate.fOptimized, corpus);
         originalPointTime += before.fPointMilliseconds;
         optimizedPointTime += after.fPointMilliseconds;
         originalRayTime += before.fRayMilliseconds;
         optimizedRayTime += after.fRayMilliseconds;
         auto &bin = FindShapeTimingBin(timingBins, candidate.fLeaves);
         ++bin.fShapeCount;
         bin.fOriginalPointMilliseconds += before.fPointMilliseconds;
         bin.fOptimizedPointMilliseconds += after.fPointMilliseconds;
         bin.fOriginalRayMilliseconds += before.fRayMilliseconds;
         bin.fOptimizedRayMilliseconds += after.fRayMilliseconds;
         std::cout << std::fixed << std::setprecision(3) << "  Timing [" << index << "] " << candidate.fName
                   << ": points " << before.fPointMilliseconds << " -> " << after.fPointMilliseconds << " ms ("
                   << (after.fPointMilliseconds > 0. ? before.fPointMilliseconds / after.fPointMilliseconds : 0.)
                   << "x), rays " << before.fRayMilliseconds << " -> " << after.fRayMilliseconds << " ms ("
                   << (after.fRayMilliseconds > 0. ? before.fRayMilliseconds / after.fRayMilliseconds : 0.)
                   << "x)\n";
      }
   }

   if (timingEnabled) {
      PrintShapeTimingBins(timingBins, pointsPerShape, raysPerShape);
      std::cout << std::fixed << std::setprecision(3)
                << "\nTGeoMultiUnion imported geometry aggregate timing (shapes=" << candidates.size()
                << ", points/shape=" << pointsPerShape << ", rays/shape=" << raysPerShape << ")\n"
                << "  Original shapes:  point containment+safety " << originalPointTime << " ms, ray distances "
                << originalRayTime << " ms\n"
                << "  Optimized shapes: point containment+safety " << optimizedPointTime << " ms, ray distances "
                << optimizedRayTime << " ms\n"
                << "  Speedup:          point containment+safety "
                << (optimizedPointTime > 0. ? originalPointTime / optimizedPointTime : 0.) << "x, ray distances "
                << (optimizedRayTime > 0. ? originalRayTime / optimizedRayTime : 0.) << "x\n";
   }
}

void ResetBooleanSelections(TGeoShape &shape)
{
   auto *composite = dynamic_cast<TGeoCompositeShape *>(&shape);
   if (!composite)
      return;
   auto *node = composite->GetBoolNode();
   node->SetSelected(0);
   ResetBooleanSelections(*node->GetLeftShape());
   ResetBooleanSelections(*node->GetRightShape());
}

NormalResponse CaptureNormal(TGeoShape &shape, const Vec3 &point, const Vec3 &direction)
{
   NormalResponse response;
   // Boolean navigation remembers the branch selected by the preceding distance query.
   // Clear that transient state so this direct normal query depends only on its arguments.
   ResetBooleanSelections(shape);
   shape.ComputeNormal(point.data(), direction.data(), response.fNormal.data());
   const Double_t length = std::sqrt(response.fNormal[0] * response.fNormal[0] +
                                     response.fNormal[1] * response.fNormal[1] +
                                     response.fNormal[2] * response.fNormal[2]);
   EXPECT_NEAR(length, 1., 1.e-12);
   return response;
}

void CompareNormal(const char *label, TGeoShape &shape, const Vec3 &point, const Vec3 &direction,
                   const NormalResponse &original)
{
   SCOPED_TRACE(label);
   const auto optimized = CaptureNormal(shape, point, direction);
   const Double_t dot = original.fNormal[0] * optimized.fNormal[0] +
                        original.fNormal[1] * optimized.fNormal[1] +
                        original.fNormal[2] * optimized.fNormal[2];
   EXPECT_GT(dot, 1. - 1.e-12) << "original normal=(" << original.fNormal[0] << ", " << original.fNormal[1]
                               << ", " << original.fNormal[2] << ") optimized normal=(" << optimized.fNormal[0]
                               << ", " << optimized.fNormal[1] << ", " << optimized.fNormal[2] << ")";
}

void ExpectLeafOnly(const TGeoMultiUnion &multiUnion)
{
   for (Int_t index = 0; index < multiUnion.GetNumberOfSolids(); ++index)
      EXPECT_FALSE(multiUnion.GetSolid(index)->IsComposite());
}

void CheckOptimizedStructure(const GeometryFixture &fixture, std::size_t depth)
{
   auto *shared = dynamic_cast<TGeoMultiUnion *>(fixture.fSharedFirst->GetShape());
   ASSERT_NE(shared, nullptr);
   EXPECT_EQ(shared->GetNumberOfSolids(), static_cast<Int_t>(depth));
   EXPECT_EQ(fixture.fSharedSecond->GetShape(), shared);
   ExpectLeafOnly(*shared);

   auto *shapeMinusUnion = dynamic_cast<TGeoCompositeShape *>(fixture.fShapeMinusUnion->GetShape());
   ASSERT_NE(shapeMinusUnion, nullptr);
   ASSERT_NE(static_cast<TGeoShape *>(shapeMinusUnion), static_cast<TGeoShape *>(fixture.fShapeMinusUnionOriginal));
   ASSERT_EQ(shapeMinusUnion->GetBoolNode()->GetBooleanOperator(), TGeoBoolNode::kGeoSubtraction);
   auto *negative = dynamic_cast<TGeoMultiUnion *>(shapeMinusUnion->GetBoolNode()->GetRightShape());
   ASSERT_NE(negative, nullptr);
   EXPECT_EQ(negative->GetNumberOfSolids(), static_cast<Int_t>(depth));
   ExpectLeafOnly(*negative);

   auto *unionMinusUnion = dynamic_cast<TGeoMultiDifference *>(fixture.fUnionMinusUnion->GetShape());
   ASSERT_NE(unionMinusUnion, nullptr);
   ASSERT_NE(static_cast<TGeoShape *>(unionMinusUnion), static_cast<TGeoShape *>(fixture.fUnionMinusUnionOriginal));
   EXPECT_EQ(unionMinusUnion->GetNpositive(), static_cast<Int_t>(depth));
   EXPECT_EQ(unionMinusUnion->GetNnegative(), static_cast<Int_t>(depth));
   ExpectLeafOnly(*unionMinusUnion);

   EXPECT_EQ(fixture.fIneligible->GetShape(), fixture.fIneligibleOriginal);
   EXPECT_EQ(fixture.fUnused->GetShape(), fixture.fUnusedOriginal);
}

} // namespace

TEST(TGeoMultiUnionGeometryEquivalence, ComplexGeometryNavigation)
{
   const TestConfig config = GetConfig();
   ASSERT_GE(config.fDepth, 3U) << "TGEOMU_DEPTH must be between 3 and 50";
   ASSERT_LE(config.fDepth, 50U) << "TGEOMU_DEPTH must be between 3 and 50";
   if (const char *filename = GeometryFile()) {
      RunImportedGeometryTest(config, filename);
      return;
   }
   TGeoManager manager("multiunion_equivalence", "multi-union geometry equivalence");
   const GeometryFixture fixture = BuildGeometry(manager, config);
   const QueryCorpus corpus = MakeQueries(config);
   PathTable paths;

   const Int_t shapeCount = manager.GetListOfShapes()->GetEntriesFast();
   ASSERT_EQ(manager.OptimizeCompositeShapes(kFALSE, 3), 3);
   ASSERT_EQ(manager.GetListOfShapes()->GetEntriesFast(), shapeCount);
   const GeometrySnapshot baseline = CaptureBaseline(manager, corpus, config, paths);
   const Bool_t timingEnabled = TimingEnabled();
   TimingResult originalTiming;
   if (timingEnabled)
      originalTiming = MeasureNavigation(manager, corpus, config);

   // These are face interiors, deliberately away from Boolean seams where more than one
   // surface normal can be valid.
   const Vec3 firstCenter = ComponentCenter(0);
   const Vec3 sharedPoint{firstCenter[0] - .5, firstCenter[1], firstCenter[2]};
   const Vec3 sharedDirection{-1., 0., 0.};
   const Vec3 subtractionPoint{7., 4., 0.};
   const Vec3 subtractionDirection{1., 0., 0.};
   const Vec3 differencePoint{firstCenter[0] - .5, firstCenter[1], firstCenter[2]};
   const Vec3 differenceDirection{-1., 0., 0.};
   const auto sharedNormal = CaptureNormal(*fixture.fSharedFirst->GetShape(), sharedPoint, sharedDirection);
   const auto subtractionNormal =
      CaptureNormal(*fixture.fShapeMinusUnion->GetShape(), subtractionPoint, subtractionDirection);
   const auto differenceNormal =
      CaptureNormal(*fixture.fUnionMinusUnion->GetShape(), differencePoint, differenceDirection);

   ASSERT_EQ(manager.OptimizeCompositeShapes(kTRUE, 3), 3);
   CheckOptimizedStructure(fixture, config.fDepth);
   ASSERT_EQ(manager.OptimizeCompositeShapes(kTRUE, 3), 0);

   CompareNormal("multi-union", *fixture.fSharedFirst->GetShape(), sharedPoint, sharedDirection, sharedNormal);
   CompareNormal("shape minus multi-union", *fixture.fShapeMinusUnion->GetShape(), subtractionPoint,
                 subtractionDirection, subtractionNormal);
   CompareNormal("multi-union minus multi-union", *fixture.fUnionMinusUnion->GetShape(), differencePoint,
                 differenceDirection, differenceNormal);

   ComparePointResponses(manager, corpus, baseline, config, paths);
   CompareRayResponses(manager, corpus, baseline, config, paths);

   if (timingEnabled) {
      const TimingResult optimizedTiming = MeasureNavigation(manager, corpus, config);
      EXPECT_EQ(optimizedTiming.fCrossings, originalTiming.fCrossings);
      PrintTimingSummary(config, originalTiming, optimizedTiming);
   }
}
