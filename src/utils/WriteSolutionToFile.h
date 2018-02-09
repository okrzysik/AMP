#ifdef USE_AMP_MESH
#ifdef USE_AMP_VECTORS
#ifndef included_WriteSolutionToFile
#define included_WriteSolutionToFile

#include "AMP/ampmesh/Mesh.h"
#include "AMP/vectors/Vector.h"

void printSolution( AMP::Mesh::Mesh::shared_ptr mesh,
                    AMP::LinearAlgebra::Vector::shared_ptr solVec,
                    std::string name );

#endif
#endif
#endif
