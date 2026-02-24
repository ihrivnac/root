// Author: Mihaela Gheata   30/03/16
/*************************************************************************
 * Copyright (C) 1995-2016, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoVGCompositeOptimiser
#define ROOT_TGeoVGCompositeOptimiser

// The class collects methods to optimise complex composite shapes using 
// MultiUnion solids where possible.
// The Recursive Dive: The function travels to the very bottom of the geometry tree.
//
// Union Consolidation: It looks for clusters of + operations. If it finds a "cloud" of 50 spheres
// joined by unions, it replaces that whole branch with a TGeoMultiUnion.

// Subtraction Efficiency: When it moves back up to a Subtraction node, the "Hole" it is subtracting
// might now be a single, fast TGeoMultiUnion instead of a messy list of 50 spheres.

#include "Rtypes.h"

#include <set>

class TGeoHMatrix;
class TGeoShape;
// class TGeoVolume;

class TGeoVGCompositeOptimiser {
private:
   // Track shapes to be deleted to avoid double-freeing shared components
   std::set<TGeoShape*> fOldShapes;
   std::set<TGeoShape*> fNewShapes;

public:
   TGeoVGCompositeOptimiser() = default;
   ~TGeoVGCompositeOptimiser() = default;

   TGeoShape* FinalRecursiveOptimize(TGeoShape* shape, Bool_t& optimised);
   // void UpdateShape(TGeoVolume* volume);
   TGeoShape* GetNewShape(TGeoShape* shape, Bool_t& optimised);

   // void VisualizeOptimization(TGeoVolume* volume);
   void Benchmark(TGeoShape* shape, Int_t nPoints = 100000);
};

#endif
