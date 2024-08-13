#ifndef included_AMP_SolveThermal
#define included_AMP_SolveThermal


#include "AMP/mesh/Mesh.h"
#include "AMP/operators/Operator.h"
#include "AMP/vectors/Vector.h"

#include <string>
#include <tuple>


namespace AMP::applications {


// Integrate a nodal or gauss point source vector to a nodal source vector
std::shared_ptr<AMP::LinearAlgebra::Vector>
integrateSourceVector( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                       std::shared_ptr<const AMP::LinearAlgebra::Vector> src,
                       std::string srcName = "",
                       std::string dstName = "" );


// Solve for the temperature
std::tuple<std::shared_ptr<AMP::LinearAlgebra::Vector>, std::shared_ptr<AMP::Operator::Operator>>
solveTemperature( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                  std::shared_ptr<const AMP::LinearAlgebra::Vector> rhs,
                  std::shared_ptr<const AMP::Database> arguments,
                  std::shared_ptr<const AMP::LinearAlgebra::Vector> initialGuess = nullptr );


// Create the non-linear coupled thermal operator
std::tuple<std::shared_ptr<AMP::Operator::Operator>, std::shared_ptr<AMP::LinearAlgebra::Vector>>
createThermalOperatorsFE( std::shared_ptr<AMP::Mesh::Mesh> mesh,
                          std::shared_ptr<const AMP::Database> input_db );


} // namespace AMP::applications

#endif
