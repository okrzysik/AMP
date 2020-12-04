

namespace AMP {
namespace LinearAlgebra {

inline ManagedMatrix::ManagedMatrix() {}

inline ManagedMatrix::ManagedMatrix( std::shared_ptr<ManagedMatrixParameters> p )
    : Matrix( p ), d_pParameters( p )
{
}

} // namespace LinearAlgebra
} // namespace AMP
