#include "AMP/vectors/data/ArrayVectorData.h"

/// \cond UNDOCUMENTED

namespace AMP::unit_test {

template<typename T>
static void testArrayVectorDimensions( std::vector<size_t> &dims, AMP::UnitTest &ut )
{
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "array" );
    auto vec = AMP::LinearAlgebra::createArrayVector<T>( dims, var );
    auto data =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::ArrayVectorData<T>>( vec->getVectorData() );
    AMP_ASSERT( data );
    auto array = data->getArray();
    if ( array.size() == ArraySize( dims ) )
        ut.passes( "ArrayVector correctly returns dimensions" );
    else
        ut.failure( "ArrayVector does not correctly returns dimensions" );
}
} // namespace AMP::unit_test
