#include "AMP/matrices/petsc/NativePetscMatrixData.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "petscmat.h"
#include "petscvec.h"


namespace AMP::LinearAlgebra {

/********************************************************
 * Constructors                                          *
 ********************************************************/
NativePetscMatrixData::NativePetscMatrixData()
{
    d_Mat                  = nullptr;
    d_MatCreatedInternally = false;
}

NativePetscMatrixData::NativePetscMatrixData( Mat m, bool internally_created )
{
    d_Mat                  = m;
    d_MatCreatedInternally = internally_created;
}

NativePetscMatrixData::NativePetscMatrixData( std::shared_ptr<MatrixParameters> params )
    : MatrixData( params )
{
    AMP_ASSERT( d_pParameters );
    const auto &comm = d_pParameters->getComm().getCommunicator();
    const auto nrows = d_pParameters->getLocalNumberOfRows();
    const auto ncols = d_pParameters->getLocalNumberOfColumns();
    MatCreate( comm, &d_Mat );
    MatSetType( d_Mat, MATMPIAIJ );
    MatSetFromOptions( d_Mat );
    MatSetSizes( d_Mat,
                 nrows,
                 ncols,
                 d_pParameters->getGlobalNumberOfRows(),
                 d_pParameters->getGlobalNumberOfColumns() );
    MatSetUp( d_Mat );
#if 1
    auto nnz_p     = d_pParameters->entryList();
    const auto nnz = *( std::max_element( nnz_p, nnz_p + nrows ) );
    MatMPIAIJSetPreallocation( d_Mat, nnz, nullptr, PETSC_DETERMINE, nullptr );
    MatSeqAIJSetPreallocation( d_Mat, nnz, nullptr );
#else
    MatMPIAIJSetPreallocation(
        d_Mat, PETSC_DEFAULT, d_pParameters->entryList(), PETSC_DEFAULT, PETSC_NULL );
    MatSeqAIJSetPreallocation( d_Mat, PETSC_DEFAULT, d_pParameters->entryList() );
#endif
    // zero out the rows explicitly

    std::vector<PetscInt> petsc_rows( nrows );
    auto rowDOFs = d_pParameters->getLeftDOFManager();
    AMP_ASSERT( rowDOFs );
    const auto srow  = rowDOFs->beginDOF();
    const auto &cols = d_pParameters->getColumns();
    std::vector<PetscInt> petsc_cols( cols.size() );
    for ( size_t i = 0; i < petsc_cols.size(); ++i ) // type conversion happening
        petsc_cols[i] = cols[i];

    auto row_nnz     = d_pParameters->entryList();
    auto current_loc = 0;

    for ( size_t i = 0; i < nrows; ++i ) {
        PetscInt global_row = srow + i;
        const auto nvals    = row_nnz[i];
        std::vector<double> vals( nvals, 0.0 );
        MatSetValues(
            d_Mat, 1, &global_row, nvals, &petsc_cols[current_loc], vals.data(), INSERT_VALUES );
        current_loc += nvals;
    }

    d_MatCreatedInternally = true;
}

NativePetscMatrixData::~NativePetscMatrixData()
{
    if ( d_MatCreatedInternally )
        PETSC::matDestroy( &d_Mat );
}


AMP_MPI NativePetscMatrixData::getComm() const
{
    if ( d_pParameters )
        return d_pParameters->getComm();

    AMP_ASSERT( d_Mat );
    MPI_Comm comm;
    PetscObjectGetComm( reinterpret_cast<PetscObject>( d_Mat ), &comm );
    return AMP_MPI( comm );
}

/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
Vector::shared_ptr NativePetscMatrixData::getRightVector() const
{
    Vec a;
    MatCreateVecs( d_Mat, &a, nullptr );
    return createVector( a, true );
}
Vector::shared_ptr NativePetscMatrixData::getLeftVector() const
{
    Vec a;
    MatCreateVecs( d_Mat, nullptr, &a );
    return createVector( a, true );
}

std::shared_ptr<Discretization::DOFManager> NativePetscMatrixData::getRightDOFManager() const
{
    if ( d_pParameters )
        return d_pParameters->getRightDOFManager();

    return getRightVector()->getDOFManager();
}

std::shared_ptr<Discretization::DOFManager> NativePetscMatrixData::getLeftDOFManager() const
{
    if ( d_pParameters )
        return d_pParameters->getLeftDOFManager();

    return getLeftVector()->getDOFManager();
}


/********************************************************
 * Get values/row by global id                           *
 ********************************************************/
size_t NativePetscMatrixData::numGlobalRows() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) rows;
}
size_t NativePetscMatrixData::numGlobalColumns() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) cols;
}
void NativePetscMatrixData::getValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) const
{
    // Zero out the data in values
    for ( size_t i = 0; i < num_rows * num_cols; i++ )
        values[i] = 0.0;
    // Get the data for each row
    auto leftDOFManager = getLeftDOFManager();
    size_t firstRow     = leftDOFManager->beginDOF();
    size_t numRows      = leftDOFManager->endDOF();

    for ( size_t i = 0; i < num_rows; i++ ) {
        if ( rows[i] < firstRow || rows[i] >= firstRow + numRows )
            continue;
        PetscInt numCols = 0;
        MatGetRow( d_Mat, rows[i], &numCols, nullptr, nullptr );
        if ( numCols == 0 )
            continue;
        MatRestoreRow( d_Mat, rows[i], &numCols, nullptr, nullptr );
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, rows[i], &numCols, &out_cols, &out_vals );
        for ( size_t j1 = 0; j1 < num_cols; j1++ ) {
            for ( int j2 = 0; j2 < numCols; j2++ ) {
                if ( (int) cols[j1] == out_cols[j2] )
                    values[i * num_cols + j1] = out_vals[j2];
            }
        }
        MatRestoreRow( d_Mat, rows[i], &numCols, &out_cols, &out_vals );
    }
}
void NativePetscMatrixData::getRowByGlobalID( size_t row,
                                              std::vector<size_t> &cols,
                                              std::vector<double> &values ) const
{
    PetscInt numCols;
    MatGetRow( d_Mat, row, &numCols, nullptr, nullptr );
    cols.resize( numCols );
    values.resize( numCols );
    MatRestoreRow( d_Mat, row, &numCols, nullptr, nullptr );
    if ( cols.size() ) { // the restore zeros out nCols
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, row, &numCols, &out_cols, &out_vals );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
        std::copy( (double *) out_vals, (double *) ( out_vals + numCols ), values.begin() );
        MatRestoreRow( d_Mat, row, &numCols, &out_cols, &out_vals );
    }
}
std::vector<size_t> NativePetscMatrixData::getColumnIDs( size_t row ) const
{
    PetscInt numCols;
    MatGetRow( d_Mat, row, &numCols, nullptr, nullptr );
    std::vector<size_t> cols( numCols );
    MatRestoreRow( d_Mat, row, &numCols, nullptr, nullptr );

    if ( cols.size() ) { // the restore zeros out nCols
        const PetscInt *out_cols;
        MatGetRow( d_Mat, row, &numCols, &out_cols, nullptr );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
        MatRestoreRow( d_Mat, row, &numCols, &out_cols, nullptr );
    }

    return cols;
}
void NativePetscMatrixData::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<PetscInt> petsc_rows( num_rows );
    std::vector<PetscInt> petsc_cols( num_cols );
    std::copy( rows, rows + num_rows, petsc_rows.begin() );
    std::copy( cols, cols + num_cols, petsc_cols.begin() );

    MatSetValues( d_Mat, num_rows, &petsc_rows[0], num_cols, &petsc_cols[0], values, ADD_VALUES );
}
void NativePetscMatrixData::setValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<PetscInt> petsc_rows( num_rows );
    std::vector<PetscInt> petsc_cols( num_cols );
    std::copy( rows, rows + num_rows, petsc_rows.begin() );
    std::copy( cols, cols + num_cols, petsc_cols.begin() );

    MatSetValues(
        d_Mat, num_rows, &petsc_rows[0], num_cols, &petsc_cols[0], values, INSERT_VALUES );
}


/********************************************************
 * Clone the matrix                                      *
 ********************************************************/
std::shared_ptr<MatrixData> NativePetscMatrixData::cloneMatrixData() const
{
    Mat new_mat;
    MatDuplicate( d_Mat, MAT_DO_NOT_COPY_VALUES, &new_mat );
    AMP_ERROR( "not quite implemented" );
    return std::make_shared<NativePetscMatrixData>( new_mat, true );
}
std::shared_ptr<MatrixData> NativePetscMatrixData::transpose() const
{
    AMP_ERROR( "Not implemented" );
}
std::shared_ptr<MatrixData> NativePetscMatrixData::duplicateMat( Mat m )
{
    Mat newMat;
    MatDuplicate( m, MAT_DO_NOT_COPY_VALUES, &newMat );
    return std::make_shared<NativePetscMatrixData>( newMat, true );
}

void NativePetscMatrixData::extractDiagonal( std::shared_ptr<Vector> v ) const
{
    auto data = std::dynamic_pointer_cast<NativePetscVectorData>( v->getVectorData() );
    MatGetDiagonal( d_Mat, data->getVec() );
}


/********************************************************
 * Copy                                                  *
 ********************************************************/
void NativePetscMatrixData::copyFromMat( Mat m ) { MatCopy( m, d_Mat, SAME_NONZERO_PATTERN ); }

/********************************************************
 * makeConsistent                                        *
 ********************************************************/
void NativePetscMatrixData::makeConsistent()
{
    MatAssemblyBegin( d_Mat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( d_Mat, MAT_FINAL_ASSEMBLY );
}


} // namespace AMP::LinearAlgebra
