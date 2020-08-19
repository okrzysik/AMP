
/// \cond UNDOCUMENTED
#include "AMP/vectors/ArrayVectorData.h"

namespace AMP {
namespace unit_test {

template<typename T>
static void testArrayVectorDimensions( std::vector<size_t> &dims, AMP::UnitTest &ut )
{
    auto var      = std::make_shared<AMP::LinearAlgebra::Variable>( "array" );
    auto vec      = AMP::LinearAlgebra::ArrayVector<T>::create( dims, var );
    auto arrayVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::ArrayVector<T>>( vec );
    AMP_ASSERT( arrayVec.get() != NULL );
    auto array = dynamic_cast<AMP::LinearAlgebra::ArrayVectorData<T>*>(arrayVec->getVectorData())->getArray();
    if ( ( array.size() == dims ) &&
         ( array.ndim() == (int) dims.size() ) )
        ut.passes( "ArrayVector correctly returns dimensions" );
    else
        ut.failure( "ArrayVector does not correctly returns dimensions" );
}
} // namespace unit_test
} // namespace AMP
