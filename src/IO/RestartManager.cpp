#include "AMP/IO/RestartManager.h"
#include "AMP/utils/UtilityMacros.h"

#include <set>


namespace AMP {


/********************************************************
 *  Write the data                                       *
 ********************************************************/
void RestartManager::write( const std::string &filename, Compression compress )
{
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    if ( globalComm.getRank() == 0 ) {
        auto fid = openHDF5( filename, "w", compress );
        closeHDF5( fid );
    }
    std::set<std::string> names;
    for ( auto &it : d_data )
        names.insert( it.first );
    globalComm.setGather( names );
    for ( auto name : names ) {
        auto it = d_data.find( name );
        int N   = globalComm.sumReduce<int>( it == d_data.end() ? 0 : 1 );
        if ( it == d_data.end() )
            continue;
        auto comm = it->second->getComm();
        if ( comm.getSize() != N )
            AMP_ERROR( "Object comm does not match registered ranks" );
        hid_t fid = 0;
        if ( comm.getRank() == 0 )
            fid = openHDF5( filename, "rw", compress );
        it->second->write( fid, filename );
        if ( comm.getRank() == 0 )
            closeHDF5( fid );
    }
}


/********************************************************
 *  Open the data                                        *
 ********************************************************/
hid_t RestartManager::open( const std::string &file ) { return openHDF5( file, "r" ); }


} // namespace AMP
