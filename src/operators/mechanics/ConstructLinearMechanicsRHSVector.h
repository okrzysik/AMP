#ifndef included_AMP_ConstructLinearMechanicsRHSVector
#define included_AMP_ConstructLinearMechanicsRHSVector

#include "AMP/materials/Material.h"

#include "AMP/utils/shared_ptr.h"

#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "AMP/ampmesh/Mesh.h"

void computeTemperatureRhsVector(
    AMP::Mesh::Mesh::shared_ptr mesh,
    AMP::shared_ptr<AMP::Database> input_db,
    AMP::LinearAlgebra::Variable::shared_ptr temperatureVar,
    AMP::LinearAlgebra::Variable::shared_ptr displacementVar,
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &currTemperatureVec,
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &prevTemperatureVec,
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec );


#endif
