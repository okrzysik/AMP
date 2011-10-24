#ifdef USE_AMP_MESH
#ifndef included_WriteSolutionToFile
#define included_WriteSolutionToFile

#include "ampmesh/MeshManager.h"

void printSolution(AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr solVec, std::string exeName);

#endif
#endif
