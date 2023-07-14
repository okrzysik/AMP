#ifndef included_TpetraVectorData_HPP_
#define included_TpetraVectorData_HPP_

#include "AMP/vectors/trilinos/tpetra/TpetraVectorData.h"


namespace AMP::LinearAlgebra {


template<typename ST, typename LO, typename GO, typename NT>
TpetraVectorData<ST, LO, GO, NT>::TpetraVectorData(
    std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
    : d_pDOFManager( dofManager )
{
    AMP_DEBUG_ASSERT( dofManager );
    const auto &mpiComm = dofManager->getComm().getCommunicator();
    auto comm           = Teuchos::rcp( new Teuchos::MpiComm<int>( mpiComm ) );
    auto map            = Teuchos::rcp( new Tpetra::Map<LO, GO, NT>(
        dofManager->numGlobalDOF(), dofManager->numLocalDOF(), comm ) );
    d_pTpetraVector     = Teuchos::rcp( new Tpetra::Vector<ST, LO, GO, NT>( map, 1 ) );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::setValuesByLocalID( size_t,
                                                           const size_t *,
                                                           const void *,
                                                           const typeID & )
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::addValuesByLocalID( size_t,
                                                           const size_t *,
                                                           const void *,
                                                           const typeID & )
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::getValuesByLocalID( size_t,
                                                           const size_t *,
                                                           void *,
                                                           const typeID & ) const
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::putRawData( const void *, const typeID & )
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::getRawData( void *, const typeID & ) const
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void *TpetraVectorData<ST, LO, GO, NT>::getRawDataBlockAsVoid( size_t )
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
const void *TpetraVectorData<ST, LO, GO, NT>::getRawDataBlockAsVoid( size_t ) const
{
    AMP_ERROR( "Not implemented" );
}

template<typename ST, typename LO, typename GO, typename NT>
void TpetraVectorData<ST, LO, GO, NT>::swapData( VectorData &other )
{
    auto otherData = dynamic_cast<TpetraVectorData<ST, LO, GO, NT> *>( &other );
    AMP_INSIST( otherData, "Not TpetraVectorData" );
    this->getTpetraVector()->swap( *( otherData->getTpetraVector() ) );
}

template<typename ST, typename LO, typename GO, typename NT>
std::shared_ptr<VectorData> TpetraVectorData<ST, LO, GO, NT>::cloneData() const
{
    return std::make_shared<TpetraVectorData<ST, LO, GO, NT>>( d_pDOFManager );
}


} // namespace AMP::LinearAlgebra
#endif
