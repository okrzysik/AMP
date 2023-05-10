#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"
#include "AMP/matrices/trilinos/EpetraMatrixData.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/data/VectorDataDefault.h"
#include "AMP/vectors/trilinos/epetra/EpetraVector.h"

#include "ProfilerApp.h"
#include <algorithm>

#include <EpetraExt_MatrixMatrix.h>
#include <Epetra_FECrsMatrix.h>
DISABLE_WARNINGS
#include <EpetraExt_Transpose_RowMatrix.h>
ENABLE_WARNINGS

#ifdef AMP_USE_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP::LinearAlgebra {

static void VerifyEpetraReturn( int err, const char *func )
{
    std::stringstream error;
    error << func << ": " << err;
    if ( err < 0 )
        AMP_ERROR( error.str() );
    if ( err > 0 )
        AMP_ERROR( error.str() );
}

/********************************************************
 * Constructors                                          *
 ********************************************************/
ManagedEpetraMatrix::ManagedEpetraMatrix( std::shared_ptr<MatrixParameters> params )
{
    d_matrixData = std::make_shared<EpetraMatrixData>( params );
}

ManagedEpetraMatrix::ManagedEpetraMatrix( const ManagedEpetraMatrix &rhs ) : Matrix( rhs )
{
    const auto rhsData = std::dynamic_pointer_cast<const EpetraMatrixData>( rhs.getMatrixData() );
    AMP_ASSERT( rhsData );
    d_matrixData = std::make_shared<EpetraMatrixData>( *rhsData );
}

ManagedEpetraMatrix::ManagedEpetraMatrix( Epetra_CrsMatrix *m, bool dele )
{
    d_matrixData = std::make_shared<EpetraMatrixData>( m, dele );
}

ManagedEpetraMatrix::ManagedEpetraMatrix( std::shared_ptr<MatrixData> data )
{
    d_matrixData = data;
}

std::shared_ptr<Matrix> ManagedEpetraMatrix::clone() const
{
    const auto data = std::dynamic_pointer_cast<const EpetraMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    auto cloneData = data->cloneMatrixData();
    return std::make_shared<ManagedEpetraMatrix>( cloneData );
}

Epetra_CrsMatrix &ManagedEpetraMatrix::getEpetra_CrsMatrix()
{
    auto data = std::dynamic_pointer_cast<EpetraMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    return data->getEpetra_CrsMatrix();
}

const Epetra_CrsMatrix &ManagedEpetraMatrix::getEpetra_CrsMatrix() const
{
    const auto data = std::dynamic_pointer_cast<const EpetraMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    return data->getEpetra_CrsMatrix();
}

std::shared_ptr<Matrix> ManagedEpetraMatrix::transpose() const
{
    EpetraExt::RowMatrix_Transpose transposer;
    Epetra_CrsMatrix &matrix = const_cast<Epetra_CrsMatrix &>( getEpetra_CrsMatrix() );
    return std::shared_ptr<Matrix>( new ManagedEpetraMatrix(
        dynamic_cast<Epetra_CrsMatrix *>( &transposer( matrix ) ), true ) );
}

/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
std::shared_ptr<Vector> ManagedEpetraMatrix::getRightVector() const
{
    const auto data = std::dynamic_pointer_cast<const EpetraMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    return data->getRightVector();
}
std::shared_ptr<Vector> ManagedEpetraMatrix::getLeftVector() const
{
    const auto data = std::dynamic_pointer_cast<const EpetraMatrixData>( d_matrixData );
    AMP_ASSERT( data );
    return data->getLeftVector();
}


/********************************************************
 * Multiply two matricies                                *
 ********************************************************/
void ManagedEpetraMatrix::multiply( shared_ptr other_op, std::shared_ptr<Matrix> &result )
{
    if ( this->numGlobalColumns() != other_op->numGlobalRows() )
        AMP_ERROR( "Inner matrix dimensions must agree" );
    if ( !std::dynamic_pointer_cast<ManagedEpetraMatrix>( other_op ) )
        AMP_ERROR( "Incompatible matrix types" );
    AMP_ASSERT( other_op->numGlobalRows() == numGlobalColumns() );
#ifdef AMP_USE_MPI
    AMP_MPI::Comm epetraComm =
        ( dynamic_cast<const Epetra_MpiComm *>( &getEpetra_CrsMatrix().RowMap().Comm() ) )->Comm();
#else
    AMP_MPI::Comm epetraComm = AMP_COMM_SELF;
#endif
    auto leftVec  = this->getLeftVector();
    auto rightVec = other_op->getRightVector();
    auto memp     = std::make_shared<MatrixParameters>(
        leftVec->getDOFManager(), rightVec->getDOFManager(), AMP_MPI( epetraComm ) );
    memp->d_CommListLeft     = leftVec->getCommunicationList();
    memp->d_CommListRight    = rightVec->getCommunicationList();
    ManagedEpetraMatrix *res = new ManagedEpetraMatrix( memp );
    PROFILE_START( "Epetra::MatrixMultiply" );
    int ierr = EpetraExt::MatrixMatrix::Multiply(
        getEpetra_CrsMatrix(),
        false,
        std::dynamic_pointer_cast<ManagedEpetraMatrix>( other_op )->getEpetra_CrsMatrix(),
        false,
        res->getEpetra_CrsMatrix(),
        true );
    AMP_ASSERT( ierr == 0 );
    PROFILE_STOP( "Epetra::MatrixMultiply" );
    result = std::shared_ptr<Matrix>( res );
}


/********************************************************
 * Multiply the matrix by a vector                       *
 ********************************************************/
void ManagedEpetraMatrix::mult( std::shared_ptr<const Vector> in, std::shared_ptr<Vector> out )
{
    AMP_ASSERT( in->getGlobalSize() == numGlobalColumns() );
    AMP_ASSERT( out->getGlobalSize() == numGlobalRows() );
    auto in_view                = EpetraVector::constView( in );
    auto out_view               = EpetraVector::view( out );
    const Epetra_Vector &in_vec = in_view->getEpetra_Vector();
    Epetra_Vector &out_vec      = out_view->getEpetra_Vector();
    int err                     = getEpetra_CrsMatrix().Multiply( false, in_vec, out_vec );
    VerifyEpetraReturn( err, "mult" );
}
void ManagedEpetraMatrix::multTranspose( std::shared_ptr<const Vector> in,
                                         std::shared_ptr<Vector> out )
{
    AMP_ASSERT( in->getGlobalSize() == numGlobalColumns() );
    AMP_ASSERT( out->getGlobalSize() == numGlobalRows() );
    auto in_view  = EpetraVector::constView( in );
    auto out_view = EpetraVector::view( out );
    int err       = getEpetra_CrsMatrix().Multiply(
        true, in_view->getEpetra_Vector(), out_view->getEpetra_Vector() );
    VerifyEpetraReturn( err, "mult" );
}


/********************************************************
 * Scale the matrix                                      *
 ********************************************************/
void ManagedEpetraMatrix::scale( double alpha )
{
    VerifyEpetraReturn( getEpetra_CrsMatrix().Scale( alpha ), "scale" );
}
void ManagedEpetraMatrix::setScalar( double ans )
{
    VerifyEpetraReturn( getEpetra_CrsMatrix().PutScalar( ans ), "setScalar" );
}
void ManagedEpetraMatrix::zero()
{
    VerifyEpetraReturn( getEpetra_CrsMatrix().PutScalar( 0.0 ), "setScalar" );
}


/********************************************************
 * axpy                                                  *
 ********************************************************/
void ManagedEpetraMatrix::axpy( double alpha, const Matrix &rhs )
{
    EpetraExt::MatrixMatrix::Add(
        dynamic_cast<const ManagedEpetraMatrix *>( &rhs )->getEpetra_CrsMatrix(),
        false,
        alpha,
        getEpetra_CrsMatrix(),
        1.0 );
}


/********************************************************
 * norm                                                  *
 ********************************************************/
double ManagedEpetraMatrix::L1Norm() const { return getEpetra_CrsMatrix().NormOne(); }


/********************************************************
 * Get/set the diagonal                                  *
 ********************************************************/
std::shared_ptr<Vector> ManagedEpetraMatrix::extractDiagonal( std::shared_ptr<Vector> vec ) const
{
    if ( !vec )
        vec = getRightVector();
    auto view = EpetraVector::view( vec );
    VerifyEpetraReturn( getEpetra_CrsMatrix().ExtractDiagonalCopy( view->getEpetra_Vector() ),
                        "extractDiagonal" );
    return vec;
}
void ManagedEpetraMatrix::setDiagonal( std::shared_ptr<const Vector> in )
{
    auto vec = EpetraVector::constView( in );
    VerifyEpetraReturn( getEpetra_CrsMatrix().ReplaceDiagonalValues( vec->getEpetra_Vector() ),
                        "setDiagonal" );
}
void ManagedEpetraMatrix::setIdentity()
{
    zero();
    int MyFirstRow = d_matrixData->getLeftDOFManager()->beginDOF();
    int MyEndRow   = d_matrixData->getLeftDOFManager()->endDOF();
    double one     = 1.0;
    for ( int i = MyFirstRow; i != MyEndRow; i++ ) {
        VerifyEpetraReturn( getEpetra_CrsMatrix().ReplaceGlobalValues( i, 1, &one, &i ),
                            "setValuesByGlobalID" );
    }
}

} // namespace AMP::LinearAlgebra
