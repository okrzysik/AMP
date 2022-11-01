#ifndef included_AMP_RestartManager_hpp
#define included_AMP_RestartManager_hpp

#include "AMP/IO/RestartManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP {


/********************************************************
 *  DataStore                                            *
 ********************************************************/
template<class TYPE>
RestartManager::DataStoreType<TYPE>::DataStoreType( std::shared_ptr<TYPE> data ) : d_data( data )
{
}
template<class TYPE>
AMP_MPI RestartManager::DataStoreType<TYPE>::getComm() const
{
    return AMP::getComm<TYPE>( *d_data );
}


/********************************************************
 *  Register data with the manager                       *
 ********************************************************/
template<class TYPE>
void RestartManager::registerData( const std::string &name, std::shared_ptr<TYPE> data )
{
    if ( d_data.find( name ) == d_data.end() )
        AMP_ERROR( name + " previously registered with restart manager" );
    d_data[name] = std::make_shared<DataStoreType<TYPE>>( data );
}


/********************************************************
 *  Load data from the manager                           *
 ********************************************************/
template<class TYPE>
std::shared_ptr<TYPE>
RestartManager::getData( hid_t fid, const std::string &name, const AMP::AMP_MPI &comm )
{
    AMP_ERROR( "Not finished" );
}


} // namespace AMP


#endif
