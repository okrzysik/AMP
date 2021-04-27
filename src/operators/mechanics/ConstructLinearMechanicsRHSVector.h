#ifndef included_AMP_ConstructLinearMechanicsRHSVector
#define included_AMP_ConstructLinearMechanicsRHSVector


#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/Database.h"
#include "AMP/vectors/Vector.h"
#include <memory>

void computeTemperatureRhsVector( AMP::Mesh::Mesh::shared_ptr mesh,
                                  std::shared_ptr<AMP::Database> input_db,
                                  AMP::LinearAlgebra::Variable::shared_ptr temperatureVar,
                                  AMP::LinearAlgebra::Variable::shared_ptr displacementVar,
                                  std::shared_ptr<AMP::LinearAlgebra::Vector> currTemperatureVec,
                                  std::shared_ptr<AMP::LinearAlgebra::Vector> prevTemperatureVec,
                                  AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


#endif
