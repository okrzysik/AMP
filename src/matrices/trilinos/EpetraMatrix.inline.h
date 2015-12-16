namespace AMP {
namespace LinearAlgebra {

inline void EpetraMatrix::VerifyEpetraReturn( int err, const char *func ) const
{
    std::stringstream error;
    error << func << ": " << err;
    if ( err < 0 ) AMP_ERROR( error.str() );
    if ( err > 0 ) AMP_ERROR( error.str() );
}

inline EpetraMatrix::EpetraMatrix( Epetra_CrsMatrix *inMatrix, bool dele )
    : d_epetraMatrix( inMatrix ), d_DeleteMatrix( dele )
{
}

inline EpetraMatrix::~EpetraMatrix()
{
    if ( d_DeleteMatrix ) delete d_epetraMatrix;
}

inline Epetra_CrsMatrix &EpetraMatrix::getEpetra_CrsMatrix() { return *d_epetraMatrix; }

inline const Epetra_CrsMatrix &EpetraMatrix::getEpetra_CrsMatrix() const { return *d_epetraMatrix; }
}
}
