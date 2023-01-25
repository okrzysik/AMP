#include "AMP/IO/RestartManager.h"
#include "AMP/IO/RestartManager.hpp"
#include "AMP/utils/Array.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"

#include <complex>
#include <set>
#include <string>
#include <vector>


namespace AMP::IO {


/********************************************************
 *  Constructor/destructor                               *
 ********************************************************/
RestartManager::RestartManager() : d_writer( true ), d_fid( -1 ) {}
RestartManager::RestartManager( const std::string &name ) : d_writer( false ), d_fid( -1 )
{
    int rank = AMP::AMP_MPI( AMP_COMM_WORLD ).getRank();
    AMP_INSIST( d_fid == hid_t( -1 ), "User must close file before opening a new one" );
    d_data.clear();
    d_names.clear();
    auto file = name + "." + AMP::Utilities::nodeToString( rank ) + ".h5";
    d_fid     = openHDF5( file, "r" );
    std::vector<std::string> names;
    std::vector<uint64_t> ids;
    readHDF5( d_fid, "RestartDataIDs", ids );
    readHDF5( d_fid, "RestartDataNames", names );
    for ( size_t i = 0; i < names.size(); i++ )
        d_names[names[i]] = ids[i];
    readCommData( name );
}
RestartManager::~RestartManager()
{
    if ( d_fid != hid_t( -1 ) )
        closeHDF5( d_fid );
    d_fid = hid_t( -1 );
}


/********************************************************
 *  Write the data                                       *
 ********************************************************/
void RestartManager::write( const std::string &name, Compression compress )
{
    int rank  = AMP::AMP_MPI( AMP_COMM_WORLD ).getRank();
    auto file = name + "." + AMP::Utilities::nodeToString( rank ) + ".h5";
    auto fid  = openHDF5( file, "w", compress );
    std::vector<std::string> names;
    std::vector<uint64_t> ids;
    for ( const auto &[name, id] : d_names ) {
        names.push_back( name );
        ids.push_back( id );
    }
    writeHDF5( fid, "RestartDataIDs", ids );
    writeHDF5( fid, "RestartDataNames", names );
    for ( const auto [id, data] : d_data ) {
        auto name = hash2String( id );
        data->write( fid, name );
    }
    closeHDF5( fid );
    writeCommData( name, compress );
}


/********************************************************
 *  Register/load communicators                          *
 ********************************************************/
void RestartManager::registerComm( const AMP::AMP_MPI &comm )
{
    auto hash = comm.hashRanks();
    if ( d_comms.find( hash ) != d_comms.end() )
        return;
    d_comms[hash] = comm;
}
AMP_MPI RestartManager::getComm( uint64_t hash )
{
    auto it = d_comms.find( hash );
    if ( it == d_comms.end() ) {
        AMP_ERROR( "Unable to find comm: " + std::to_string( hash ) );
    }
    return it->second.dup();
}
void RestartManager::writeCommData( const std::string &name, Compression compress )
{
    // Collect the comm data
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    int rank = globalComm.getRank();
    int size = globalComm.getSize();
    std::set<uint64_t> comm_ids;
    for ( auto &[id, comm] : d_comms )
        comm_ids.insert( id );
    globalComm.setGather( comm_ids );
    std::vector<uint64_t> comm_list( comm_ids.begin(), comm_ids.end() );
    std::vector<bool> comm_data( size * comm_ids.size(), false );
    for ( size_t i = 0; i < comm_ids.size(); i++ ) {
        std::vector<bool> tmp( size, false );
        if ( d_comms.find( comm_list[i] ) != d_comms.end() )
            tmp[rank] = true;
        globalComm.anyReduce( tmp );
        for ( int j = 0; j < size; j++ ) {
            if ( tmp[j] )
                comm_data[i * size + j] = true;
        }
    }
    if ( rank == 0 ) {
        auto file = name + ".comms.h5";
        auto fid  = openHDF5( file, "w", compress );
        writeHDF5( fid, "ids", comm_list );
        writeHDF5( fid, "data", comm_data );
        closeHDF5( fid );
    }
}
void RestartManager::readCommData( const std::string &name )
{
    auto file = name + ".comms.h5";
    auto fid  = openHDF5( file, "r" );
    std::vector<uint64_t> comm_list;
    std::vector<bool> comm_data;
    readHDF5( fid, "ids", comm_list );
    readHDF5( fid, "data", comm_data );
    closeHDF5( fid );
    d_comms.clear();
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    int rank = globalComm.getRank();
    int size = globalComm.getSize();
    for ( size_t i = 0; i < comm_list.size(); i++ ) {
        bool test = comm_data[i * size + rank];
        auto comm = globalComm.split( test ? 1 : 0 );
        if ( test )
            d_comms[comm_list[i]] = comm;
    }
}


/********************************************************
 *  Register data with the manager                       *
 ********************************************************/
void RestartManager::registerData( std::shared_ptr<DataStore> data )
{
    auto hash    = data->getHash();
    d_data[hash] = data;
}


/********************************************************
 *  Register data with the manager                       *
 ********************************************************/
std::string AMP::IO::RestartManager::hash2String( uint64_t id )
{
    char id_chars[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789#$";
    std::string name;
    while ( id > 0 ) {
        name += id_chars[id & 0x3F];
        id >>= 6;
    }
    return name;
}


/********************************************************
 *  Explicit instantiations                              *
 ********************************************************/
template<class TYPE>
AMP::IO::RestartManager::DataStoreType<TYPE>::DataStoreType( const std::string &name,
                                                             std::shared_ptr<const TYPE> data,
                                                             RestartManager * )
    : d_name( name ), d_data( data )
{
}
template<class TYPE>
uint64_t AMP::IO::RestartManager::DataStoreType<TYPE>::getHash() const
{
    return AMP::Utilities::hash_char( d_name.data() );
}
template<class TYPE>
void RestartManager::DataStoreType<TYPE>::write( hid_t fid, const std::string &name ) const
{
    writeHDF5( fid, name, *d_data );
}
template<class TYPE>
std::shared_ptr<TYPE> RestartManager::getData( uint64_t hash )
{
    AMP_INSIST( d_fid != hid_t( -1 ), "User must open file before reading data" );
    auto data = std::make_shared<TYPE>();
    readHDF5( d_fid, hash2String( hash ), *data );
    return data;
}
#define INSTANTIATE( TYPE )                                                          \
    template class RestartManager::DataStoreType<TYPE>;                              \
    template void RestartManager::registerData( const std::string &, const TYPE & ); \
    template std::shared_ptr<TYPE> RestartManager::getData( const std::string & )
#define INSTANTIATE2( TYPE )          \
    INSTANTIATE( TYPE );              \
    INSTANTIATE( std::vector<TYPE> ); \
    INSTANTIATE( AMP::Array<TYPE> )
INSTANTIATE( bool );
INSTANTIATE( char );
INSTANTIATE( uint8_t );
INSTANTIATE( uint16_t );
INSTANTIATE( uint32_t );
INSTANTIATE( uint64_t );
INSTANTIATE( int8_t );
INSTANTIATE( int16_t );
INSTANTIATE( int32_t );
INSTANTIATE( int64_t );
INSTANTIATE( float );
INSTANTIATE( double );
INSTANTIATE( std::complex<float> );
INSTANTIATE( std::complex<double> );
INSTANTIATE( std::byte );
INSTANTIATE( std::string );
INSTANTIATE( std::string_view );
INSTANTIATE( AMP::Database );


} // namespace AMP::IO
