

namespace AMP {
namespace LinearAlgebra {

inline NativeVector::NativeVector() {}

inline NativeVector::NativeVector( parameters_ptr params )
    : Vector( std::dynamic_pointer_cast<VectorParameters>( params ) )
{
}
} // namespace LinearAlgebra
} // namespace AMP
