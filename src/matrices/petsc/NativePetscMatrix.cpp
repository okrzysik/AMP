#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "petscmat.h"
#include "petscvec.h"


namespace AMP {
namespace LinearAlgebra {


// Get vector
static std::shared_ptr<Vec> getVec( std::shared_ptr<Vector> v )
{
    auto data = std::dynamic_pointer_cast<NativePetscVectorData>( v->getVectorData() );
    if ( data )
        return std::shared_ptr<Vec>( new Vec( data->getVec() ), []( auto ) {} );
    return std::shared_ptr<Vec>( new Vec( PETSC::getVec( v ) ), []( Vec *v ) { VecDestroy( v ); } );
}
static std::shared_ptr<Vec> getVec( std::shared_ptr<const Vector> v )
{
    return getVec( std::const_pointer_cast<Vector>( v ) );
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
NativePetscMatrix::NativePetscMatrix( Mat m, bool internally_created )
{
    d_Mat                  = m;
    d_MatCreatedInternally = internally_created;
    MPI_Comm comm;
    PetscObjectGetComm( reinterpret_cast<PetscObject>( m ), &comm );
    d_comm = AMP_MPI( comm );
}
NativePetscMatrix::NativePetscMatrix() { d_MatCreatedInternally = false; }
NativePetscMatrix::~NativePetscMatrix()
{
    if ( d_MatCreatedInternally )
        PETSC::matDestroy( &d_Mat );
}


/********************************************************
 * Multiply the matrix by a vector                       *
 ********************************************************/
void NativePetscMatrix::multiply( shared_ptr other_op, shared_ptr &result )
{
    auto other = std::dynamic_pointer_cast<NativePetscMatrix>( other_op );
    AMP_INSIST( other != nullptr, "Incompatible matrix types" );

    auto res = new NativePetscMatrix();
    MatMatMult( d_Mat, other->d_Mat, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &res->d_Mat );
    result.reset( res );
}


/********************************************************
 * Get the diagonal                                      *
 ********************************************************/
Vector::shared_ptr NativePetscMatrix::extractDiagonal( Vector::shared_ptr v ) const
{
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<NativePetscVectorData>( v->getVectorData() ) ) {
        retVal = v;
    } else {
        retVal = getRightVector();
        retVal->setVariable( v->getVariable() );
    }
    auto data = std::dynamic_pointer_cast<NativePetscVectorData>( v->getVectorData() );
    MatGetDiagonal( d_Mat, data->getVec() );
    return retVal;
}


/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
Vector::shared_ptr NativePetscMatrix::getRightVector() const
{
    Vec a;
#if PETSC_VERSION_LE( 3, 2, 0 )
    MatGetVecs( d_Mat, &a, PETSC_NULL );
#else
    MatCreateVecs( d_Mat, &a, PETSC_NULL );
#endif
    return createVector( a, true );
}
Vector::shared_ptr NativePetscMatrix::getLeftVector() const
{
    Vec a;
#if PETSC_VERSION_LE( 3, 2, 0 )
    MatGetVecs( d_Mat, PETSC_NULL, &a );
#else
    MatCreateVecs( d_Mat, PETSC_NULL, &a );
#endif
    return createVector( a, true );
}
Discretization::DOFManager::shared_ptr NativePetscMatrix::getRightDOFManager() const
{
    return getRightVector()->getDOFManager();
}
Discretization::DOFManager::shared_ptr NativePetscMatrix::getLeftDOFManager() const
{
    return getLeftVector()->getDOFManager();
}


/********************************************************
 * Get values/row by global id                           *
 ********************************************************/
size_t NativePetscMatrix::numGlobalRows() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) rows;
}
size_t NativePetscMatrix::numGlobalColumns() const
{
    int rows, cols;
    MatGetSize( d_Mat, &rows, &cols );
    return (size_t) cols;
}
void NativePetscMatrix::getValuesByGlobalID(
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
        int numCols = 0;
        MatGetRow( d_Mat, rows[i], &numCols, PETSC_NULL, PETSC_NULL );
        if ( numCols == 0 )
            continue;
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, rows[i], &numCols, &out_cols, &out_vals );
        for ( size_t j1 = 0; j1 < num_cols; j1++ ) {
            for ( int j2 = 0; j2 < numCols; j2++ ) {
                if ( (int) cols[j1] == out_cols[j2] )
                    values[i * num_cols + j1] = out_vals[j2];
            }
        }
    }
}
void NativePetscMatrix::getRowByGlobalID( size_t row,
                                          std::vector<size_t> &cols,
                                          std::vector<double> &values ) const
{
    int numCols;
    MatGetRow( d_Mat, row, &numCols, PETSC_NULL, PETSC_NULL );
    cols.resize( numCols );
    values.resize( numCols );
    if ( numCols ) {
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, row, &numCols, &out_cols, &out_vals );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
        std::copy( (double *) out_vals, (double *) ( out_vals + numCols ), values.begin() );
    }
}
std::vector<size_t> NativePetscMatrix::getColumnIDs( size_t row ) const
{
    int numCols;
    MatGetRow( d_Mat, row, &numCols, PETSC_NULL, PETSC_NULL );
    std::vector<size_t> cols( numCols );

    if ( numCols ) {
        const PetscInt *out_cols;
        const PetscScalar *out_vals;
        MatGetRow( d_Mat, row, &numCols, &out_cols, &out_vals );
        std::copy(
            (unsigned int *) out_cols, (unsigned int *) ( out_cols + numCols ), cols.begin() );
    }

    return cols;
}
void NativePetscMatrix::addValuesByGlobalID(
    size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values )
{
    std::vector<PetscInt> petsc_rows( num_rows );
    std::vector<PetscInt> petsc_cols( num_cols );
    std::copy( rows, rows + num_rows, petsc_rows.begin() );
    std::copy( cols, cols + num_cols, petsc_cols.begin() );

    MatSetValues( d_Mat, num_rows, &petsc_rows[0], num_cols, &petsc_cols[0], values, ADD_VALUES );
}
void NativePetscMatrix::setValuesByGlobalID(
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
std::shared_ptr<Matrix> NativePetscMatrix::cloneMatrix() const
{
    Mat new_mat;
    MatDuplicate( d_Mat, MAT_DO_NOT_COPY_VALUES, &new_mat );
    AMP_ERROR( "not quite implemented" );
    return shared_ptr( new NativePetscMatrix( new_mat, true ) );
}
std::shared_ptr<Matrix> NativePetscMatrix::duplicateMat( Mat m )
{
    Mat newMat;
    MatDuplicate( m, MAT_DO_NOT_COPY_VALUES, &newMat );
    return std::shared_ptr<Matrix>( new NativePetscMatrix( newMat, true ) );
}


/********************************************************
 * Copy                                                  *
 ********************************************************/
void NativePetscMatrix::copyFromMat( Mat m ) { MatCopy( m, d_Mat, SAME_NONZERO_PATTERN ); }


/********************************************************
 * Multiply                                              *
 ********************************************************/
void NativePetscMatrix::mult( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMult( d_Mat, *getVec( in ), *getVec( out ) );
}
void NativePetscMatrix::multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMultTranspose( d_Mat, *getVec( in ), *getVec( out ) );
}


/********************************************************
 * axpy                                                  *
 ********************************************************/
void NativePetscMatrix::axpy( double alpha, const Matrix &x )
{
    AMP_ASSERT( x.numGlobalRows() == this->numGlobalRows() );
    AMP_ASSERT( x.numGlobalColumns() == this->numGlobalColumns() );
    MatAXPY(
        d_Mat, alpha, dynamic_cast<const NativePetscMatrix *>( &x )->d_Mat, SAME_NONZERO_PATTERN );
}


/********************************************************
 * Scale                                                 *
 ********************************************************/
void NativePetscMatrix::scale( double alpha ) { MatScale( d_Mat, alpha ); }


/********************************************************
 * Set to scalar                                         *
 ********************************************************/
void NativePetscMatrix::setScalar( double ans )
{
    if ( ans != 0.0 )
        AMP_ERROR( "Cannot perform operation on NativePetscMatrix yet!" );
    MatZeroEntries( d_Mat );
}
void NativePetscMatrix::zero() { MatZeroEntries( d_Mat ); }


/********************************************************
 * setDiagonal                                           *
 ********************************************************/
void NativePetscMatrix::setDiagonal( Vector::const_shared_ptr in )
{
    MatDiagonalSet( d_Mat, *getVec( in ), INSERT_VALUES );
}
void NativePetscMatrix::setIdentity()
{
    MatZeroEntries( d_Mat );
    MatShift( d_Mat, 1.0 );
}

/********************************************************
 * makeConsistent                                        *
 ********************************************************/
void NativePetscMatrix::makeConsistent()
{
    MatAssemblyBegin( d_Mat, MAT_FINAL_ASSEMBLY );
    MatAssemblyEnd( d_Mat, MAT_FINAL_ASSEMBLY );
}


/********************************************************
 * Norm                                                  *
 ********************************************************/
double NativePetscMatrix::L1Norm() const
{
    double retVal;
    MatNorm( d_Mat, NORM_1, &retVal );
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
