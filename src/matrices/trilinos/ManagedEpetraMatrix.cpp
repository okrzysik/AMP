#include "AMP/matrices/trilinos/ManagedEpetraMatrix.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/trilinos/epetra/EpetraVector.h"

#include "ProfilerApp.h"
#include <algorithm>

#include <EpetraExt_MatrixMatrix.h>
#include <Epetra_FECrsMatrix.h>
DISABLE_WARNINGS
#include <EpetraExt_Transpose_RowMatrix.h>
ENABLE_WARNINGS

#ifdef USE_EXT_MPI
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

namespace AMP::LinearAlgebra {


static inline auto createEpetraMap( std::shared_ptr<AMP::Discretization::DOFManager> DOFs,
                                    const AMP_MPI &comm )
{
#ifdef USE_EXT_MPI
    Epetra_MpiComm comm2 = comm.getCommunicator();
#else
    NULL_USE( comm );
    Epetra_SerialComm comm2;
#endif
    AMP_INSIST( DOFs->numGlobalDOF() < 0x80000000,
                "Epetra does not support vectors with global size greater than 2^31" );
    auto N_local  = static_cast<int>( DOFs->numLocalDOF() );
    auto N_global = static_cast<int>( DOFs->numGlobalDOF() );
    return std::make_shared<Epetra_Map>( N_global, N_local, 0, comm2 );
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
ManagedEpetraMatrix::ManagedEpetraMatrix( std::shared_ptr<ManagedMatrixParameters> params )
    : EpetraMatrix( *createEpetraMap( params->getLeftDOFManager(), params->getComm() ),
                    nullptr,
                    params->entryList() ),
      ManagedMatrix( params )
{
    AMP_ASSERT( !d_comm.isNull() );
}
ManagedEpetraMatrix::ManagedEpetraMatrix( const ManagedEpetraMatrix &rhs )
    : Matrix(),
      EpetraMatrix(
          *createEpetraMap( rhs.d_pParameters->getLeftDOFManager(), rhs.d_pParameters->getComm() ),
          nullptr,
          rhs.d_pParameters->entryList() ),
      ManagedMatrix( rhs.d_pParameters )
{
    for ( size_t i = d_pParameters->getLeftDOFManager()->beginDOF();
          i != d_pParameters->getLeftDOFManager()->endDOF();
          i++ ) {
        std::vector<size_t> cols;
        std::vector<double> vals;
        rhs.getRowByGlobalID( (int) i, cols, vals );
        std::vector<size_t> cols2( cols.size() );
        for ( size_t j = 0; j != cols.size(); j++ )
            cols2[j] = cols[j];
        createValuesByGlobalID( i, cols2 );
    }
    d_RangeMap  = rhs.d_RangeMap;
    d_DomainMap = rhs.d_DomainMap;
    makeConsistent();
}
ManagedEpetraMatrix::ManagedEpetraMatrix( Epetra_CrsMatrix *m, bool dele )
    : EpetraMatrix( m, dele ), ManagedMatrix( nullptr )
{
}
std::shared_ptr<Matrix> ManagedEpetraMatrix::cloneMatrix() const
{
    auto *r           = new ManagedEpetraMatrix( *this );
    r->d_DeleteMatrix = true;
    return std::shared_ptr<Matrix>( r );
}


/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
std::shared_ptr<Vector> ManagedEpetraMatrix::getRightVector() const
{
    int localSize  = d_pParameters->getLocalNumberOfColumns();
    int globalSize = d_pParameters->getGlobalNumberOfColumns();
    int localStart = d_pParameters->getRightDOFManager()->beginDOF();
    auto buffer    = std::make_shared<VectorDataCPU<double>>( localStart, localSize, globalSize );
    auto vec       = createEpetraVector(
        d_pParameters->d_CommListRight, d_pParameters->getRightDOFManager(), buffer );
    vec->setVariable( d_pParameters->d_VariableRight );
    return vec;
}
std::shared_ptr<Vector> ManagedEpetraMatrix::getLeftVector() const
{
    int localSize  = d_pParameters->getLocalNumberOfRows();
    int globalSize = d_pParameters->getGlobalNumberOfRows();
    int localStart = d_pParameters->getRightDOFManager()->beginDOF();
    auto buffer    = std::make_shared<VectorDataCPU<double>>( localStart, localSize, globalSize );
    auto vec       = createEpetraVector(
        d_pParameters->d_CommListLeft, d_pParameters->getLeftDOFManager(), buffer );
    vec->setVariable( d_pParameters->d_VariableLeft );
    return vec;
}
Discretization::DOFManager::shared_ptr ManagedEpetraMatrix::getRightDOFManager() const
{
    return d_pParameters->getRightDOFManager();
}
Discretization::DOFManager::shared_ptr ManagedEpetraMatrix::getLeftDOFManager() const
{
    return d_pParameters->getLeftDOFManager();
}
size_t ManagedEpetraMatrix::numGlobalRows() const { return d_epetraMatrix->NumGlobalRows(); }
size_t ManagedEpetraMatrix::numGlobalColumns() const { return d_epetraMatrix->NumGlobalCols(); }


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
#ifdef USE_EXT_MPI
    MPI_Comm epetraComm =
        ( dynamic_cast<const Epetra_MpiComm *>( &d_epetraMatrix->RowMap().Comm() ) )->Comm();
#else
    MPI_Comm epetraComm = AMP_COMM_SELF;
#endif
    auto leftVec  = this->getLeftVector();
    auto rightVec = other_op->getRightVector();
    auto memp     = std::make_shared<ManagedMatrixParameters>(
        leftVec->getDOFManager(), rightVec->getDOFManager(), AMP_MPI( epetraComm ) );
    memp->d_CommListLeft     = leftVec->getCommunicationList();
    memp->d_CommListRight    = rightVec->getCommunicationList();
    ManagedEpetraMatrix *res = new ManagedEpetraMatrix( memp );
    PROFILE_START( "Epetra::MatrixMultiply" );
    int ierr = EpetraExt::MatrixMatrix::Multiply(
        *d_epetraMatrix,
        false,
        *( std::dynamic_pointer_cast<ManagedEpetraMatrix>( other_op )->d_epetraMatrix ),
        false,
        *( res->d_epetraMatrix ),
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
    int err                     = d_epetraMatrix->Multiply( false, in_vec, out_vec );
    VerifyEpetraReturn( err, "mult" );
}
void ManagedEpetraMatrix::multTranspose( std::shared_ptr<const Vector> in,
                                         std::shared_ptr<Vector> out )
{
    AMP_ASSERT( in->getGlobalSize() == numGlobalColumns() );
    AMP_ASSERT( out->getGlobalSize() == numGlobalRows() );
    auto in_view  = EpetraVector::constView( in );
    auto out_view = EpetraVector::view( out );
    int err =
        d_epetraMatrix->Multiply( true, in_view->getEpetra_Vector(), out_view->getEpetra_Vector() );
    VerifyEpetraReturn( err, "mult" );
}


/********************************************************
 * Scale the matrix                                      *
 ********************************************************/
void ManagedEpetraMatrix::scale( double alpha )
{
    VerifyEpetraReturn( d_epetraMatrix->Scale( alpha ), "scale" );
}
void ManagedEpetraMatrix::setScalar( double ans )
{
    VerifyEpetraReturn( d_epetraMatrix->PutScalar( ans ), "setScalar" );
}
void ManagedEpetraMatrix::zero()
{
    VerifyEpetraReturn( d_epetraMatrix->PutScalar( 0.0 ), "setScalar" );
}


/********************************************************
 * axpy                                                  *
 ********************************************************/
void ManagedEpetraMatrix::axpy( double alpha, const Matrix &rhs )
{
    EpetraExt::MatrixMatrix::Add(
        *( dynamic_cast<const ManagedEpetraMatrix *>( &rhs )->d_epetraMatrix ),
        false,
        alpha,
        *d_epetraMatrix,
        1.0 );
}


/********************************************************
 * norm                                                  *
 ********************************************************/
double ManagedEpetraMatrix::L1Norm() const { return d_epetraMatrix->NormOne(); }


/********************************************************
 * setOtherData                                          *
 ********************************************************/
void ManagedEpetraMatrix::setOtherData()
{
    AMP_MPI myComm = d_pParameters->getComm();
    int ndxLen     = d_OtherData.size();
    int totNdxLen  = myComm.sumReduce( ndxLen );
    if ( totNdxLen == 0 ) {
        return;
    }
    int dataLen  = 0;
    auto cur_row = d_OtherData.begin();
    while ( cur_row != d_OtherData.end() ) {
        dataLen += cur_row->second.size();
        ++cur_row;
    }
    auto rows   = new int[dataLen + 1]; // Add one to have the new work
    auto cols   = new int[dataLen + 1];
    auto data   = new double[dataLen + 1];
    int cur_ptr = 0;
    cur_row     = d_OtherData.begin();
    while ( cur_row != d_OtherData.end() ) {
        auto cur_elem = cur_row->second.begin();
        while ( cur_elem != cur_row->second.end() ) {
            rows[cur_ptr] = cur_row->first;
            cols[cur_ptr] = cur_elem->first;
            data[cur_ptr] = cur_elem->second;
            ++cur_ptr;
            ++cur_elem;
        }
        ++cur_row;
    }

    int totDataLen = myComm.sumReduce( dataLen );

    auto aggregateRows = new int[totDataLen];
    auto aggregateCols = new int[totDataLen];
    auto aggregateData = new double[totDataLen];

    myComm.allGather( rows, dataLen, aggregateRows );
    myComm.allGather( cols, dataLen, aggregateCols );
    myComm.allGather( data, dataLen, aggregateData );

    int MyFirstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    int MyEndRow   = d_pParameters->getLeftDOFManager()->endDOF();
    for ( int i = 0; i != totDataLen; i++ ) {
        if ( ( aggregateRows[i] >= MyFirstRow ) && ( aggregateRows[i] < MyEndRow ) ) {
            setValueByGlobalID( aggregateRows[i], aggregateCols[i], aggregateData[i] );
        }
    }

    d_OtherData.clear();
    delete[] rows;
    delete[] cols;
    delete[] data;
    delete[] aggregateRows;
    delete[] aggregateCols;
    delete[] aggregateData;
}


void ManagedEpetraMatrix::fillComplete() { EpetraMatrix::fillComplete(); }


/********************************************************
 * Get/set the diagonal                                  *
 ********************************************************/
std::shared_ptr<Vector> ManagedEpetraMatrix::extractDiagonal( std::shared_ptr<Vector> vec ) const
{
    if ( !vec )
        vec = getRightVector();
    auto view = EpetraVector::view( vec );
    VerifyEpetraReturn( d_epetraMatrix->ExtractDiagonalCopy( view->getEpetra_Vector() ),
                        "extractDiagonal" );
    return vec;
}
void ManagedEpetraMatrix::setDiagonal( std::shared_ptr<const Vector> in )
{
    auto vec = EpetraVector::constView( in );
    VerifyEpetraReturn( d_epetraMatrix->ReplaceDiagonalValues( vec->getEpetra_Vector() ),
                        "setDiagonal" );
}
void ManagedEpetraMatrix::setIdentity()
{
    zero();
    int MyFirstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    int MyEndRow   = d_pParameters->getLeftDOFManager()->endDOF();
    double one     = 1.0;
    for ( int i = MyFirstRow; i != MyEndRow; i++ ) {
        VerifyEpetraReturn( d_epetraMatrix->ReplaceGlobalValues( i, 1, &one, &i ),
                            "setValuesByGlobalID" );
    }
}


/********************************************************
 * Set/Add values by global id                           *
 ********************************************************/
void ManagedEpetraMatrix::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<int> epetra_cols( num_cols );
    std::copy( cols, cols + num_cols, epetra_cols.begin() );

    for ( size_t i = 0; i != num_rows; i++ )
        VerifyEpetraReturn( d_epetraMatrix->SumIntoGlobalValues(
                                rows[i], num_cols, values + num_cols * i, epetra_cols.data() ),
                            "addValuesByGlobalId" );
}
void ManagedEpetraMatrix::createValuesByGlobalID( size_t row, const std::vector<size_t> &cols )
{
    if ( cols.empty() )
        return;
    std::vector<int> indices( cols.size() );
    std::copy( cols.begin(), cols.end(), indices.begin() );

    std::vector<double> values( cols.size(), 0 );

    VerifyEpetraReturn(
        d_epetraMatrix->InsertGlobalValues( row, cols.size(), values.data(), indices.data() ),
        "setValuesByGlobalID" );
}
void ManagedEpetraMatrix::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<int> epetra_cols( num_cols );
    std::copy( cols, cols + num_cols, epetra_cols.begin() );

    size_t MyFirstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    size_t MyEndRow   = d_pParameters->getLeftDOFManager()->endDOF();
    for ( size_t i = 0; i != num_rows; i++ ) {
        VerifyEpetraReturn( d_epetraMatrix->ReplaceGlobalValues(
                                rows[i], num_cols, values + num_cols * i, epetra_cols.data() ),
                            "setValuesByGlobalID" );
        if ( rows[i] < MyFirstRow || rows[i] >= MyEndRow ) {
            for ( size_t j = 0; j != num_cols; j++ ) {
                d_OtherData[rows[i]][cols[j]] = values[num_cols * i + j];
            }
        }
    }
}


/********************************************************
 * Get values/row by global id                           *
 ********************************************************/
void ManagedEpetraMatrix::getValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const
{
    // Zero out the data in values
    for ( size_t i = 0; i < num_rows * num_cols; i++ )
        values[i] = 0.0;
    // Get the data for each row
    size_t firstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    size_t numRows  = d_pParameters->getLeftDOFManager()->endDOF();
    std::vector<int> row_cols;
    std::vector<double> row_values;
    for ( size_t i = 0; i < num_rows; i++ ) {
        if ( rows[i] < firstRow || rows[i] >= firstRow + numRows )
            continue;
        size_t localRow = rows[i] - firstRow;
        int numCols     = d_pParameters->entriesInRow( localRow );
        if ( numCols == 0 )
            continue;
        row_cols.resize( numCols );
        row_values.resize( numCols );
        VerifyEpetraReturn( d_epetraMatrix->ExtractGlobalRowCopy(
                                rows[i], numCols, numCols, &( row_values[0] ), &( row_cols[0] ) ),
                            "getValuesByGlobalID" );
        for ( size_t j1 = 0; j1 < num_cols; j1++ ) {
            for ( size_t j2 = 0; j2 < (size_t) numCols; j2++ ) {
                if ( cols[j1] == (size_t) row_cols[j2] )
                    values[i * num_cols + j1] = row_values[j2];
            }
        }
    }
}
void ManagedEpetraMatrix::getRowByGlobalID( size_t row,
                                            std::vector<size_t> &cols,
                                            std::vector<double> &values ) const
{
    size_t firstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    size_t numRows  = d_pParameters->getLeftDOFManager()->endDOF();
    AMP_ASSERT( row >= firstRow );
    AMP_ASSERT( row < firstRow + numRows );

    size_t localRow = row - firstRow;
    int numCols     = d_pParameters->entriesInRow( localRow );
    cols.resize( numCols );
    values.resize( numCols );

    if ( numCols ) {
        std::vector<int> epetra_cols( numCols );
        VerifyEpetraReturn( d_epetraMatrix->ExtractGlobalRowCopy(
                                row, numCols, numCols, &( values[0] ), &( epetra_cols[0] ) ),
                            "getRowByGlobalID" );
        std::copy( epetra_cols.begin(), epetra_cols.end(), cols.begin() );
    }
}


std::vector<size_t> ManagedEpetraMatrix::getColumnIDs( size_t row ) const
{
    size_t firstRow = d_pParameters->getLeftDOFManager()->beginDOF();
    size_t numRows  = d_pParameters->getLeftDOFManager()->endDOF();
    AMP_ASSERT( row >= firstRow );
    AMP_ASSERT( row < firstRow + numRows );

    size_t localRow = row - firstRow;
    int numCols     = d_pParameters->entriesInRow( localRow );
    std::vector<size_t> cols( numCols );

    if ( numCols ) {
        std::vector<double> values( numCols );
        std::vector<int> epetra_cols( numCols );
        VerifyEpetraReturn( d_epetraMatrix->ExtractGlobalRowCopy(
                                row, numCols, numCols, &( values[0] ), &( epetra_cols[0] ) ),
                            "getRowByGlobalID" );
        std::copy( epetra_cols.begin(), epetra_cols.end(), cols.begin() );
    }

    return cols;
}

/********************************************************
 * makeConsistent                                        *
 ********************************************************/
void ManagedEpetraMatrix::makeConsistent()
{
    auto *mat = dynamic_cast<Epetra_FECrsMatrix *>( d_epetraMatrix );
    if ( mat ) {
        VerifyEpetraReturn( mat->GlobalAssemble( false ), "makeParallelConsistent" );
        fillComplete();
    }
    setOtherData();
}


void ManagedEpetraMatrix::FillComplete() { d_epetraMatrix->FillComplete(); }

} // namespace AMP::LinearAlgebra
