#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/matrices/petsc/NativePetscMatrix.h"
#include "AMP/vectors/ExternalVectorDeleter.h"
#include "AMP/vectors/petsc/NativePetscVector.h"


namespace AMP {
namespace Operator {

PetscMatrixShellOperator::PetscMatrixShellOperator(
    const AMP::shared_ptr<OperatorParameters> &params )
    : LinearOperator( params )
{
    d_iMatLocalRowSize    = -1234;
    d_iMatLocalColumnSize = -5678;
}

void PetscMatrixShellOperator::setComm( AMP_MPI comm ) { d_comm = comm; }

void PetscMatrixShellOperator::setMatLocalRowSize( int val ) { d_iMatLocalRowSize = val; }

void PetscMatrixShellOperator::setMatLocalColumnSize( int val ) { d_iMatLocalColumnSize = val; }

void PetscMatrixShellOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                      AMP::LinearAlgebra::Vector::shared_ptr f )
{
    d_operator->apply( u, f );
}

void PetscMatrixShellOperator::reset( const AMP::shared_ptr<OperatorParameters> &params )
{
    d_operator->reset( params );
}

AMP::LinearAlgebra::Variable::shared_ptr PetscMatrixShellOperator::getOutputVariable()
{
    return d_operator->getOutputVariable();
}

AMP::LinearAlgebra::Variable::shared_ptr PetscMatrixShellOperator::getInputVariable()
{
    return d_operator->getInputVariable();
}

void PetscMatrixShellOperator::setOperator( AMP::shared_ptr<Operator> op )
{
    d_operator = op;
    MatCreateShell( d_comm.getCommunicator(),
                    d_iMatLocalRowSize,
                    d_iMatLocalColumnSize,
                    PETSC_DETERMINE,
                    PETSC_DETERMINE,
                    this,
                    &d_mat );
    MatShellSetOperation(
        d_mat, MATOP_MULT, ( void ( * )() )( &( PetscMatrixShellOperator::mult ) ) );
    d_matrix.reset( new AMP::LinearAlgebra::NativePetscMatrix( d_mat, true ) );
}

PetscErrorCode PetscMatrixShellOperator::mult( Mat mat, Vec in, Vec out )
{
    void *ctx;
    MatShellGetContext( mat, &ctx );
    auto *op = reinterpret_cast<PetscMatrixShellOperator *>( ctx );

    AMP::LinearAlgebra::Vector::shared_ptr inVec(
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( in->data ),
        AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::LinearAlgebra::Vector::shared_ptr outVec(
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( out->data ),
        AMP::LinearAlgebra::ExternalVectorDeleter() );

    ( op->d_operator )->apply( inVec, outVec );

    return ( 0 );
}
} // namespace Operator
} // namespace AMP
