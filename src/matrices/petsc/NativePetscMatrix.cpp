#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/petsc/NativePetscMatrixData.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "petscmat.h"
#include "petscvec.h"


namespace AMP::LinearAlgebra {


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

static Mat getMat( Matrix *m )
{
    AMP_ASSERT( m );
    auto data = dynamic_cast<NativePetscMatrixData *>( m->getMatrixData().get() );
    AMP_ASSERT( data );
    return data->getMat();
}
static Mat getMat( const Matrix *m ) { return getMat( const_cast<const Matrix *>( m ) ); }

static Mat getMat( std::shared_ptr<const Matrix> m )
{
    return getMat( std::const_pointer_cast<Matrix>( m ) );
}

static Mat getMat( std::shared_ptr<Matrix> m )
{
    auto data = std::dynamic_pointer_cast<NativePetscMatrixData>( m->getMatrixData() );
    AMP_ASSERT( data );
    return data->getMat();
}


/********************************************************
 * Constructors                                          *
 ********************************************************/
NativePetscMatrix::NativePetscMatrix( Mat m, bool internally_created )
{
    d_matrixData = std::make_shared<NativePetscMatrixData>( m, internally_created );
}
NativePetscMatrix::~NativePetscMatrix() {}


/********************************************************
 * Multiply the matrix by a vector                       *
 ********************************************************/
void NativePetscMatrix::multiply( shared_ptr other_op, shared_ptr &result )
{
    auto other = std::dynamic_pointer_cast<NativePetscMatrix>( other_op );
    AMP_INSIST( other != nullptr, "Incompatible matrix types" );

    Mat resMat;
    MatMatMult( getMat( this ), getMat( other_op ), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &resMat );

    auto res = new NativePetscMatrix( resMat, true ); // see if it should be false...
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
    MatGetDiagonal( getMat( this ), data->getVec() );
    return retVal;
}


/********************************************************
 * Get the left/right Vector/DOFManager                  *
 ********************************************************/
Vector::shared_ptr NativePetscMatrix::getRightVector() const
{
    std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->getRightVector();
}
Vector::shared_ptr NativePetscMatrix::getLeftVector() const
{
    std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->getLeftVector();
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
 * Clone the matrix                                      *
 ********************************************************/
std::shared_ptr<Matrix> NativePetscMatrix::cloneMatrix() const
{
    Mat new_mat;
    MatDuplicate( getMat( this ), MAT_DO_NOT_COPY_VALUES, &new_mat );
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
void NativePetscMatrix::copyFromMat( Mat m )
{
    std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->copyFromMat( m );
}


/********************************************************
 * Multiply                                              *
 ********************************************************/
void NativePetscMatrix::mult( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMult( getMat( this ), *getVec( in ), *getVec( out ) );
}
void NativePetscMatrix::multTranspose( Vector::const_shared_ptr in, Vector::shared_ptr out )
{
    MatMultTranspose( getMat( this ), *getVec( in ), *getVec( out ) );
}


/********************************************************
 * axpy                                                  *
 ********************************************************/
void NativePetscMatrix::axpy( double alpha, const Matrix &x )
{
    AMP_ASSERT( x.numGlobalRows() == this->numGlobalRows() );
    AMP_ASSERT( x.numGlobalColumns() == this->numGlobalColumns() );
    MatAXPY( getMat( this ), alpha, getMat( &x ), SAME_NONZERO_PATTERN );
}


/********************************************************
 * Scale                                                 *
 ********************************************************/
void NativePetscMatrix::scale( double alpha ) { MatScale( getMat( this ), alpha ); }


/********************************************************
 * Set to scalar                                         *
 ********************************************************/
void NativePetscMatrix::setScalar( double ans )
{
    if ( ans != 0.0 )
        AMP_ERROR( "Cannot perform operation on NativePetscMatrix yet!" );
    MatZeroEntries( getMat( this ) );
}
void NativePetscMatrix::zero() { MatZeroEntries( getMat( this ) ); }


/********************************************************
 * setDiagonal                                           *
 ********************************************************/
void NativePetscMatrix::setDiagonal( Vector::const_shared_ptr in )
{
    MatDiagonalSet( getMat( this ), *getVec( in ), INSERT_VALUES );
}
void NativePetscMatrix::setIdentity()
{
    auto mat = getMat( this );
    MatZeroEntries( mat );
    MatShift( mat, 1.0 );
}

/********************************************************
 * Norm                                                  *
 ********************************************************/
double NativePetscMatrix::L1Norm() const
{
    double retVal;
    MatNorm( getMat( this ), NORM_1, &retVal );
    return retVal;
}


} // namespace AMP::LinearAlgebra
