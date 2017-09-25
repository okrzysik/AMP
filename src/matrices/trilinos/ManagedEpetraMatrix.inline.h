namespace AMP {
namespace LinearAlgebra {


inline ManagedEpetraMatrix::ManagedEpetraMatrix( Epetra_CrsMatrix *m, bool dele )
    : EpetraMatrix( m, dele ), ManagedMatrix( MatrixParameters::shared_ptr() )
{
}


inline double ManagedEpetraMatrix::L1Norm() const { return d_epetraMatrix->NormOne(); }


inline Matrix::shared_ptr ManagedEpetraMatrix::cloneMatrix() const
{
    ManagedEpetraMatrix *r = new ManagedEpetraMatrix( *this );
    r->d_DeleteMatrix      = true;
    return Matrix::shared_ptr( r );
}
} // namespace LinearAlgebra
} // namespace AMP
