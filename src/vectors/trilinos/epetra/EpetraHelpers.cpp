#include "AMP/vectors/trilinos/epetra/EpetraHelpers.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"


namespace AMP::LinearAlgebra {


static inline double *getBufferPtr( std::shared_ptr<VectorData> buf )
{
    size_t N_blocks = buf->numberOfDataBlocks();
    if ( N_blocks == 0 )
        return nullptr;
    if ( N_blocks > 1 )
        AMP_ERROR( "More than 1 data block detected" );
    return buf->getRawDataBlock<double>( 0 );
}


/********************************************************
 * Get an Epetra vector from an AMP vector               *
 ********************************************************/
std::shared_ptr<Epetra_Vector> getEpetra( std::shared_ptr<Vector> vec )
{
    EpetraVectorEngineParameters params(
        vec->getLocalSize(), vec->getComm(), vec->getCommunicationList() );
    auto ptr  = getBufferPtr( vec->getVectorData() );
    auto vec2 = std::make_shared<Epetra_Vector>( View, params.getEpetraMap(), ptr );
    return vec2;
}


} // namespace AMP::LinearAlgebra
