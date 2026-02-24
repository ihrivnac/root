
// How the Logic Flows
// The Recursive Dive: The function travels to the very bottom of the geometry tree.

// Union Consolidation: It looks for clusters of + operations. If it finds a "cloud" of 50 spheres
// joined by unions, it replaces that whole branch with a TGeoMultiUnion.

// Subtraction Efficiency: When it moves back up to a Subtraction node, the "Hole" it is subtracting
// might now be a single, fast TGeoMultiUnion instead of a messy list of 50 spheres.

// Usage:

// TGeoVGCompositeOptimiser opt;
// TGeoVolume* myVol = gGeoManager->GetVolume("MY_VOLUME_NAME");

// // This is the "Best Results" entry point
// opt.UpdateShape(myVol); 

// // Verification
// opt.VisualizeOptimization(myVol);
// opt.Benchmark(myVol->GetShape());

#include "TGeoVGCompositeOptimiser.h"

#include "TGeoBoolNode.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoVGMultiUnion.h"
#include "TStopwatch.h"

#include <set>
#include <vector>

namespace {

////////////////////////////////////////////////////////////////////////////////
/// Helper to check if all components of composite shape are VecGeom solids

void CheckMembers(TGeoShape* shape, Bool_t& isVecGeom)
{    
    if (!shape->IsComposite()) {
        return;
    }

    TGeoCompositeShape* comp = (TGeoCompositeShape*)shape;
    TGeoBoolNode* node = comp->GetBoolNode();

    TGeoShape* leftShape = node->GetLeftShape();
    if ((! leftShape->IsComposite()) && (!leftShape->IsVecGeom())) {
        isVecGeom = false;
        std::cerr << "Component shape " << leftShape->GetName() << " is not VecGeom" << std::endl;
        return;
    } 
    CheckMembers(leftShape, isVecGeom);

    TGeoShape* rightShape = node->GetRightShape();
    if ((! rightShape->IsComposite()) && (! rightShape->IsVecGeom())) {
        isVecGeom = false;
        std::cerr << "Component shape " << rightShape->GetName() << " is not VecGeom" << std::endl;
        return;
    } 
}

////////////////////////////////////////////////////////////////////////////////
/// Helper to collect components for a MultiUnion

void CollectMembers(TGeoShape* shape, TGeoHMatrix* currentMat, 
       std::vector<TGeoShape*>& shapes, std::vector<TGeoHMatrix*>& matrices, 
       TGeoBoolNode::EGeoBoolType targetType)
{    
    if (!shape->IsComposite()) {
        shapes.push_back(shape);
        matrices.push_back(new TGeoHMatrix(*currentMat));
        return;
    }

    TGeoCompositeShape* comp = (TGeoCompositeShape*)shape;
    TGeoBoolNode* node = comp->GetBoolNode();

    if (node->GetBooleanOperator() == targetType) {
        // Recurse Left
        TGeoHMatrix leftMat = (*currentMat) * (*node->GetLeftMatrix());
        CollectMembers(node->GetLeftShape(), &leftMat, shapes, matrices, targetType);
        
        // Recurse Right
        TGeoHMatrix rightMat = (*currentMat) * (*node->GetRightMatrix());
        CollectMembers(node->GetRightShape(), &rightMat, shapes, matrices, targetType);
    } else {
        // We hit a different operation (e.g., hit a Subtraction while looking for Unions)
        // We add this optimized branch as a single component.
        shapes.push_back(shape);
        matrices.push_back(new TGeoHMatrix(*currentMat));
    }
}

}

////////////////////////////////////////////////////////////////////////////////
/// Run optimisation

TGeoShape* TGeoVGCompositeOptimiser::FinalRecursiveOptimize(TGeoShape* shape,  Bool_t& optimised) 
{
    if (!shape || !shape->IsComposite()) return shape;

    // Record this shape as one to be replaced/cleaned
    fOldShapes.insert(shape);

    TGeoCompositeShape* comp = (TGeoCompositeShape*)shape;
    TGeoBoolNode* node = comp->GetBoolNode();
    
    // 1. Process children first
    TGeoShape* leftOpt  = FinalRecursiveOptimize(node->GetLeftShape(), optimised);
    TGeoShape* rightOpt = FinalRecursiveOptimize(node->GetRightShape(), optimised);

    TGeoShape* result = nullptr;

    // 2. Apply MultiUnion Flattening for Unions
    if (node->GetBooleanOperator() == TGeoBoolNode::kGeoUnion) {
        std::vector<TGeoShape*> shapes;
        std::vector<TGeoHMatrix*> matrices;
        TGeoHMatrix identity;
        CollectMembers(shape, &identity, shapes, matrices, TGeoBoolNode::kGeoUnion);

        if (shapes.size() > 2) {
            std::vector<ShapeEntry> components;
            for (size_t i = 0; i < shapes.size(); ++i) {
                components.push_back({shapes[i], matrices[i]});
            }
            auto multi = new TGeoVGMultiUnion(Form("%s_MU", shape->GetName()), components);
            // for (size_t i = 0; i < shapes.size(); ++i) multi->AddNode(shapes[i], matrices[i]);
            // multi->Voxelize();
            fNewShapes.insert(multi);
            optimised = true;
            result = multi;
        }
    }

    // 3. Fallback: Create a new Composite with the already optimized children
    if (!result) {
        TGeoBoolNode* newNode = nullptr;
        if (node->GetBooleanOperator() == TGeoBoolNode::kGeoUnion)
            newNode = new TGeoUnion(leftOpt, rightOpt, (TGeoMatrix*)node->GetLeftMatrix(), (TGeoMatrix*)node->GetRightMatrix());
        else if (node->GetBooleanOperator() == TGeoBoolNode::kGeoSubtraction)
            newNode = new TGeoSubtraction(leftOpt, rightOpt, (TGeoMatrix*)node->GetLeftMatrix(), (TGeoMatrix*)node->GetRightMatrix());
        else
            newNode = new TGeoIntersection(leftOpt, rightOpt, (TGeoMatrix*)node->GetLeftMatrix(), (TGeoMatrix*)node->GetRightMatrix());
        
        result = new TGeoCompositeShape(Form("%s_opt", shape->GetName()), newNode);
        fNewShapes.insert(result);
    }

    gGeoManager->AddShape(result);
    return result;
}

// void TGeoVGCompositeOptimiser::UpdateShape(TGeoVolume* volume) {
//     if (!volume) return;
    
//     TGeoShape* oldShape = volume->GetShape();
//     if (!oldShape->IsComposite()) return;

//     // Clear the tracking set
//     fOldShapes.clear();

//     // Generate optimized shape
//     TGeoShape* newShape = FinalRecursiveOptimize(oldShape);

//     // Update the volume to use the new shape
//     volume->SetShape(newShape);

//     // Remove old shapes from TGeoManager to free memory
//     for (auto* s : fOldShapes) {
//         gGeoManager->GetListOfShapes()->Remove(s);
//         // Note: Actual 'delete s' can be risky if shapes are shared across volumes.
//         // It is often safer to let TGeoManager handle the list and just Remove().
//     }
    
//     std::cout << "Optimization complete. Old shapes removed from Manager." << std::endl;
// }

TGeoShape* TGeoVGCompositeOptimiser::GetNewShape(TGeoShape* shape, Bool_t& optimised) 
{
    if (!shape) {
        optimised = false;
        return nullptr;
    }

    TGeoShape* oldShape = shape;
    if (!oldShape->IsComposite()) {
        optimised = false;
        return oldShape;
    }

    Bool_t isVecGeom = true;
    CheckMembers(shape, isVecGeom);
    if ( ! isVecGeom) {
        optimised = false;
        return oldShape;                
    }

    // Clear the tracking set
    fOldShapes.clear();
    fNewShapes.clear();

    // Generate optimized shape
    TGeoShape* newShape = FinalRecursiveOptimize(oldShape, optimised);

    if (optimised) {
       std::cout << "Shape was optimised." << std::endl;
       // Remove old shapes from TGeoManager to free memory
       for (auto* s : fOldShapes) {
           gGeoManager->GetListOfShapes()->Remove(s);
           // Note: Actual 'delete s' can be risky if shapes are shared across volumes.
           // It is often safer to let TGeoManager handle the list and just Remove().
       }
      std::cout << "Optimization complete. Old shapes removed from Manager." << std::endl;
      return newShape;

    }
    else {
       std::cout << "Shape was not optimised." << std::endl;
       // Remove new shapes from TGeoManager to free memory
       for (auto* s : fNewShapes) {
           gGeoManager->GetListOfShapes()->Remove(s);
           // Note: Actual 'delete s' can be risky if shapes are shared across volumes.
           // It is often safer to let TGeoManager handle the list and just Remove().
       }
      std::cout << "Optimization not performed. New shapes removed from Manager." << std::endl;
      return shape;
    }
}

// void TGeoVGCompositeOptimiser::VisualizeOptimization(TGeoVolume* vol) {
//     // 1. Draw the volume in a GL Viewer
//     vol->Draw("ogl");

//     // 2. To see the Voxelization (The "Magic" behind MultiUnion speed)
//     TGeoShape* shape = vol->GetShape();
//     if (shape->IsA() == TGeoMultiUnion::Class()) {
//         TGeoMultiUnion* mu = (TGeoMultiUnion*)shape;
        
//         // This helper displays the internal voxel structure in the terminal
//         mu->InspectShape(); 
        
//         // To see voxels in the 3D view, we use the GeoManager painter
//         gGeoManager->GetGeomPainter()->SetCheckMode(kTRUE);
//         vol->Draw("voxels");
//     }
// }

void TGeoVGCompositeOptimiser::Benchmark(TGeoShape* shape, Int_t nPoints) {
    Double_t point[3] = {0, 0, 0};
    Double_t dir[3] = {1, 0, 0};
    
    TStopwatch timer;
    timer.Start();
    for (int i=0; i<nPoints; ++i) {
        shape->DistFromOutside(point, dir);
    }
    timer.Stop();
    
    printf("Time for %d checks: %f seconds\n", nPoints, timer.RealTime());
}

