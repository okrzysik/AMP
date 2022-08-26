#ifndef included_AMP_WriteSolutionToFile
#define included_AMP_WriteSolutionToFile

#include "AMP/mesh/Mesh.h"
#include "AMP/vectors/Vector.h"

void printSolution( AMP::Mesh::Mesh::shared_ptr mesh,
                    AMP::LinearAlgebra::Vector::shared_ptr solVec,
                    const std::string &name );

#endif
