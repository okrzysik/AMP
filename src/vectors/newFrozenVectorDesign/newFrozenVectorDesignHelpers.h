
#ifndef included_AMP_NewFrozenVectorDesignHelpers
#define included_AMP_NewFrozenVectorDesignHelpers

#include "vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {

AMP::LinearAlgebra::Vector::shared_ptr
subsetExceptForVariable( AMP::LinearAlgebra::Vector::shared_ptr inVec,
                         AMP::LinearAlgebra::Variable::shared_ptr var );

AMP::LinearAlgebra::Vector::shared_ptr joinVectors( AMP::LinearAlgebra::Vector::shared_ptr vec1,
                                                    AMP::LinearAlgebra::Vector::shared_ptr vec2 );
}
}

#endif
