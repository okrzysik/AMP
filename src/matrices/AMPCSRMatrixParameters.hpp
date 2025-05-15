#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/AMPCSRMatrixParameters.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"

namespace AMP::LinearAlgebra {

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    AMP::Utilities::Backend backend,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, backend ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<Variable> varLeft,
    std::shared_ptr<Variable> varRight,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, varLeft, varRight ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<Variable> varLeft,
    std::shared_ptr<Variable> varRight,
    AMP::Utilities::Backend backend,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, varLeft, varRight, backend ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<CommunicationList> commListLeft,
    std::shared_ptr<CommunicationList> commListRight,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, commListLeft, commListRight ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<CommunicationList> commListLeft,
    std::shared_ptr<CommunicationList> commListRight,
    AMP::Utilities::Backend backend,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, commListLeft, commListRight, backend ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<Variable> varLeft,
    std::shared_ptr<Variable> varRight,
    std::shared_ptr<CommunicationList> commListLeft,
    std::shared_ptr<CommunicationList> commListRight,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, varLeft, varRight, commListLeft, commListRight ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

template<typename Policy>
AMPCSRMatrixParameters<Policy>::AMPCSRMatrixParameters(
    std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
    std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
    const AMP_MPI &comm,
    std::shared_ptr<Variable> varLeft,
    std::shared_ptr<Variable> varRight,
    std::shared_ptr<CommunicationList> commListLeft,
    std::shared_ptr<CommunicationList> commListRight,
    AMP::Utilities::Backend backend,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::lidx_t &,
                              typename Policy::lidx_t & )> getRowNNZ,
    const std::function<void( const typename Policy::gidx_t,
                              typename Policy::gidx_t *,
                              typename Policy::gidx_t * )> getRowCols )
    : MatrixParameters( dofLeft, dofRight, comm, varLeft, varRight, commListLeft, commListRight, backend ),
      d_getRowNNZ( getRowNNZ ),
      d_getRowCols( getRowCols )
{
    if ( d_getRowNNZ || d_getRowCols ) {
        AMP_INSIST(
            d_getRowNNZ && d_getRowCols,
            "AMPCSRMatrixParameters: Must either provide both getRow functions or neither" );
    }
}

} // namespace AMP::LinearAlgebra
