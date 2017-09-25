
/// \cond UNDOCUMENTED

namespace AMP {
namespace unit_test {

template<typename T>
void testArrayVectorDimensions( std::vector<size_t> &dims, AMP::UnitTest &ut )
{
    AMP::LinearAlgebra::Variable::shared_ptr var( new AMP::LinearAlgebra::Variable( "array" ) );
    AMP::LinearAlgebra::Vector::shared_ptr vec =
        AMP::LinearAlgebra::ArrayVector<T>::create( dims, var );
    AMP::shared_ptr<AMP::LinearAlgebra::ArrayVector<T>> arrayVec =
        std::dynamic_pointer_cast<AMP::LinearAlgebra::ArrayVector<T>>( vec );
    AMP_ASSERT( arrayVec.get() != NULL );
    if ( ( arrayVec->getArray().size() == dims ) &&
         ( arrayVec->getArray().ndim() == (int) dims.size() ) )
        ut.passes( "ArrayVector correctly returns dimensions" );
    else
        ut.failure( "ArrayVector does not correctly returns dimensions" );
}
} // namespace unit_test
} // namespace AMP
