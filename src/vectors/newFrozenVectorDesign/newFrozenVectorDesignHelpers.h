
#ifndef included_AMP_NewFrozenVectorDesignHelpers
#define included_AMP_NewFrozenVectorDesignHelpers

#include "AMP/vectors/Vector.h"

namespace AMP::LinearAlgebra {

AMP::LinearAlgebra::Vector::shared_ptr
subsetExceptForVariable( AMP::LinearAlgebra::Vector::shared_ptr inVec,
                         std::shared_ptr<AMP::LinearAlgebra::Variable> var );

AMP::LinearAlgebra::Vector::shared_ptr joinVectors( AMP::LinearAlgebra::Vector::shared_ptr vec1,
                                                    AMP::LinearAlgebra::Vector::shared_ptr vec2 );
} // namespace AMP::LinearAlgebra

#endif
