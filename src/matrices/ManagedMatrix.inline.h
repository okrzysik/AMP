

namespace AMP::LinearAlgebra {

inline ManagedMatrix::ManagedMatrix() {}

inline ManagedMatrix::ManagedMatrix( std::shared_ptr<ManagedMatrixParameters> p )
    : Matrix( p ), d_pParameters( p )
{
    d_comm = p->getComm();
    AMP_ASSERT( !d_comm.isNull() );
}

} // namespace AMP::LinearAlgebra
