#include "utils/AMP_MPI.h"
#include "vectors/ArrayVector.h"

/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {

template <typename T>
bool
testArrayVectorDimensions( std::vector<size_t> &dims)
{
    AMP::LinearAlgebra::Variable::shared_ptr var(new AMP::LinearAlgebra::Variable( "array" ));
    AMP::LinearAlgebra::Vector::shared_ptr vec = AMP::LinearAlgebra::ArrayVector<T>::create(dims, var );
    auto arrayVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::ArrayVector<T>>(vec);
    return ((arrayVec->getArray().size()==dims)&&(arrayVec->getArray().ndim()==dims.size()));
}

}
}
