#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/petsc/NativePetscMatrixData.h"
#include "AMP/matrices/petsc/NativePetscMatrixOperations.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/petsc/NativePetscVectorData.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "petscmat.h"

namespace AMP::LinearAlgebra {

NativePetscMatrix::NativePetscMatrix()
{
    d_matrixData = std::make_shared<NativePetscMatrixData>();
    d_matrixOps  = std::make_shared<NativePetscMatrixOperations>();
}

NativePetscMatrix::NativePetscMatrix( Mat m, bool internally_created )
{
    d_matrixData = std::make_shared<NativePetscMatrixData>( m, internally_created );
    d_matrixOps  = std::make_shared<NativePetscMatrixOperations>();
}

NativePetscMatrix::NativePetscMatrix( std::shared_ptr<MatrixParametersBase> params,
				      const std::function<std::vector<size_t>( size_t )> getRow  )
    : Matrix( params )
{
  d_matrixData = std::make_shared<NativePetscMatrixData>( params, getRow );
    d_matrixOps  = std::make_shared<NativePetscMatrixOperations>();
}

NativePetscMatrix::NativePetscMatrix( std::shared_ptr<MatrixData> data ) : Matrix( data )
{
    d_matrixOps = std::make_shared<NativePetscMatrixOperations>();
}

NativePetscMatrix::~NativePetscMatrix() {}

void NativePetscMatrix::multiply( shared_ptr other_op, shared_ptr &result )
{
    auto other = std::dynamic_pointer_cast<NativePetscMatrix>( other_op );
    AMP_INSIST( other != nullptr, "Incompatible matrix types" );

    result = std::make_shared<NativePetscMatrix>();

    d_matrixOps->matMultiply(
        *d_matrixData, *( other_op->getMatrixData() ), *result->getMatrixData() );
}

Vector::shared_ptr NativePetscMatrix::extractDiagonal( Vector::shared_ptr v ) const
{
    Vector::shared_ptr retVal = v;
    if ( !retVal ) {
        retVal = this->getRightVector();
    } else if ( std::dynamic_pointer_cast<NativePetscVectorData>( v->getVectorData() ) ) {
        retVal = v;
    } else {
        AMP_ERROR( "NativePetscMatrix::extractDiagonal(): Not handled for vectors that are not "
                   "NativePetscVector" );
    }

    d_matrixData->extractDiagonal( retVal );
    return retVal;
}

Vector::shared_ptr NativePetscMatrix::getRightVector() const
{
    return std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->getRightVector();
}
Vector::shared_ptr NativePetscMatrix::getLeftVector() const
{
    return std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->getLeftVector();
}
std::shared_ptr<Discretization::DOFManager> NativePetscMatrix::getRightDOFManager() const
{
    //    return getRightVector()->getDOFManager();
    return d_matrixData->getRightDOFManager();
}
std::shared_ptr<Discretization::DOFManager> NativePetscMatrix::getLeftDOFManager() const
{
    return d_matrixData->getLeftDOFManager();
    //    return getLeftVector()->getDOFManager();
}

std::shared_ptr<Matrix> NativePetscMatrix::clone() const
{
    return std::make_shared<NativePetscMatrix>( d_matrixData->cloneMatrixData() );
}

std::shared_ptr<Matrix> NativePetscMatrix::duplicateMat( Mat m )
{
    auto data = NativePetscMatrixData::duplicateMat( m );
    return std::make_shared<NativePetscMatrix>( data );
}

void NativePetscMatrix::copyFromMat( Mat m )
{
    std::dynamic_pointer_cast<NativePetscMatrixData>( d_matrixData )->copyFromMat( m );
}

} // namespace AMP::LinearAlgebra
