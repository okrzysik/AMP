#ifndef included_AMP_WriteSolutionToFile
#define included_AMP_WriteSolutionToFile

#include "AMP/mesh/Mesh.h"
#include "AMP/vectors/Vector.h"

void printSolution( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                    AMP::LinearAlgebra::Vector::shared_ptr solVec,
                    const std::string &name );

#endif
