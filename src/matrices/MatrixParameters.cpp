#include "AMP/matrices/MatrixParameters.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/CommunicationList.h"
#include "AMP/vectors/Vector.h"


namespace AMP::LinearAlgebra {


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );

    // comm lists not provided, generate from dof managers
    generateCommLists();
}

MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    AMP::Utilities::Backend backend,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, backend ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );

    // comm lists not provided, generate from dof managers
    generateCommLists();
}

MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<Variable> varLeft,
                                    std::shared_ptr<Variable> varRight,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, varLeft, varRight ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );

    // comm lists not provided, generate from dof managers
    generateCommLists();
}


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<Variable> varLeft,
                                    std::shared_ptr<Variable> varRight,
                                    AMP::Utilities::Backend backend,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, varLeft, varRight, backend ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );

    // comm lists not provided, generate from dof managers
    generateCommLists();
}

MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<CommunicationList> commListLeft,
                                    std::shared_ptr<CommunicationList> commListRight,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_CommListLeft( commListLeft ),
      d_CommListRight( commListRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );
}


MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<CommunicationList> commListLeft,
                                    std::shared_ptr<CommunicationList> commListRight,
                                    AMP::Utilities::Backend backend,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, backend ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_CommListLeft( commListLeft ),
      d_CommListRight( commListRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );
}

MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<Variable> varLeft,
                                    std::shared_ptr<Variable> varRight,
                                    std::shared_ptr<CommunicationList> commListLeft,
                                    std::shared_ptr<CommunicationList> commListRight,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, varLeft, varRight ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_CommListLeft( commListLeft ),
      d_CommListRight( commListRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );
}

MatrixParameters::MatrixParameters( std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
                                    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
                                    const AMP_MPI &comm,
                                    std::shared_ptr<Variable> varLeft,
                                    std::shared_ptr<Variable> varRight,
                                    std::shared_ptr<CommunicationList> commListLeft,
                                    std::shared_ptr<CommunicationList> commListRight,
                                    AMP::Utilities::Backend backend,
                                    const std::function<std::vector<size_t>( size_t )> getRow )
    : MatrixParametersBase( comm, varLeft, varRight, backend ),
      d_DOFManagerLeft( dofLeft ),
      d_DOFManagerRight( dofRight ),
      d_CommListLeft( commListLeft ),
      d_CommListRight( commListRight ),
      d_getRow( getRow )
{
    AMP_ASSERT( d_DOFManagerLeft );
    AMP_ASSERT( d_DOFManagerRight );
}

void MatrixParameters::generateCommLists()
{
    // parameters objects filled in from dof managers
    auto leftClParams            = std::make_shared<CommunicationListParameters>();
    leftClParams->d_comm         = d_comm;
    leftClParams->d_localsize    = d_DOFManagerLeft->numLocalDOF();
    leftClParams->d_remote_DOFs  = d_DOFManagerLeft->getRemoteDOFs();
    auto rightClParams           = std::make_shared<CommunicationListParameters>();
    rightClParams->d_comm        = d_comm;
    rightClParams->d_localsize   = d_DOFManagerRight->numLocalDOF();
    rightClParams->d_remote_DOFs = d_DOFManagerRight->getRemoteDOFs();

    // comm lists from those
    d_CommListLeft  = std::make_shared<CommunicationList>( leftClParams );
    d_CommListRight = std::make_shared<CommunicationList>( rightClParams );
}

size_t MatrixParameters::getLocalNumberOfRows() const { return d_DOFManagerLeft->numLocalDOF(); }

size_t MatrixParameters::getLocalNumberOfColumns() const
{
    return d_DOFManagerRight->numLocalDOF();
}

size_t MatrixParameters::getGlobalNumberOfRows() const { return d_DOFManagerLeft->numGlobalDOF(); }

size_t MatrixParameters::getGlobalNumberOfColumns() const
{
    return d_DOFManagerRight->numGlobalDOF();
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getLeftDOFManager()
{
    return d_DOFManagerLeft;
}

std::shared_ptr<AMP::Discretization::DOFManager> MatrixParameters::getRightDOFManager()
{
    return d_DOFManagerRight;
}

std::shared_ptr<CommunicationList> MatrixParameters::getLeftCommList() { return d_CommListLeft; }

std::shared_ptr<CommunicationList> MatrixParameters::getRightCommList() { return d_CommListRight; }

} // namespace AMP::LinearAlgebra
