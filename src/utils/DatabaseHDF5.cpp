#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/utils/Database.h"
#include "AMP/utils/FactoryStrategy.hpp"


using KeyDataFactory = AMP::FactoryStrategy<AMP::KeyData>;


/****************************************************************************
 * read/write HDF5                                                           *
 ****************************************************************************/
#ifdef AMP_USE_HDF5
void AMP::Database::writeHDF5( int64_t fid, std::string_view name ) const
{
    hid_t gid = createGroup( fid, name );
    AMP::writeHDF5( gid, "check", d_check );
    AMP::writeHDF5( gid, "name", d_name );
    AMP::writeHDF5( gid, "keys", d_keys );
    std::vector<std::string> types( d_data.size() );
    for ( size_t i = 0; i < d_data.size(); i++ )
        types[i] = d_data[i]->getClassType().name;
    AMP::writeHDF5( gid, "types", types );
    for ( size_t i = 0; i < d_data.size(); i++ )
        d_data[i]->writeHDF5( gid, d_keys[i] );
    closeGroup( gid );
}
void AMP::Database::readHDF5( int64_t fid, std::string_view name )
{
    hid_t gid = openGroup( fid, name );
    std::vector<std::string> types;
    AMP::readHDF5( gid, "check", d_check );
    AMP::readHDF5( gid, "name", d_name );
    AMP::readHDF5( gid, "keys", d_keys );
    AMP::readHDF5( gid, "types", types );
    d_data.clear();
    d_data.resize( d_keys.size() );
    d_hash.resize( d_keys.size(), 0 );
    d_used.resize( d_keys.size(), false );
    for ( size_t i = 0; i < d_data.size(); i++ ) {
        d_data[i] = KeyDataFactory::create( types[i] );
        d_data[i]->readHDF5( gid, d_keys[i] );
        d_hash[i] = hashString( d_keys[i] );
    }
    closeGroup( gid );
}
#else
void AMP::Database::writeHDF5( int64_t, std::string_view ) const {}
void AMP::Database::readHDF5( int64_t, std::string_view ) {}
#endif


/************************************************************************
 * read/write HDF5 (AMP::Database)                                      *
 ***********************************************************************/
#ifdef AMP_USE_HDF5
template<>
hid_t AMP::getHDF5datatype<AMP::Database>()
{
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::writeHDF5Array<AMP::Database>( hid_t,
                                         const std::string_view &,
                                         const AMP::Array<AMP::Database> & )
{
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::readHDF5Array<AMP::Database>( hid_t,
                                        const std::string_view &,
                                        AMP::Array<AMP::Database> & )
{
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::writeHDF5Scalar<AMP::Database>( hid_t fid,
                                          const std::string_view &name,
                                          const AMP::Database &data )
{
    data.writeHDF5( fid, name );
}
template<>
void AMP::readHDF5Scalar<AMP::Database>( hid_t fid,
                                         const std::string_view &name,
                                         AMP::Database &data )
{
    data.readHDF5( fid, name );
}
#endif


/************************************************************************
 * read/write HDF5 (AMP::Database::Check)                               *
 ***********************************************************************/
#ifdef AMP_USE_HDF5
template<>
hid_t AMP::getHDF5datatype<AMP::Database::Check>()
{
    return getHDF5datatype<uint8_t>();
}
template<>
void AMP::writeHDF5Array<AMP::Database::Check>( hid_t,
                                                const std::string_view &,
                                                const AMP::Array<AMP::Database::Check> & )
{
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::readHDF5Array<AMP::Database::Check>( hid_t,
                                               const std::string_view &,
                                               AMP::Array<AMP::Database::Check> & )
{
    AMP_ERROR( "Not finished" );
}
template<>
void AMP::writeHDF5Scalar<AMP::Database::Check>( hid_t fid,
                                                 const std::string_view &name,
                                                 const AMP::Database::Check &data )
{
    auto data2 = static_cast<uint8_t>( data );
    AMP::writeHDF5( fid, name, data2 );
}
template<>
void AMP::readHDF5Scalar<AMP::Database::Check>( hid_t fid,
                                                const std::string_view &name,
                                                AMP::Database::Check &data )
{
    uint8_t data2;
    AMP::readHDF5( fid, name, data2 );
    data = static_cast<AMP::Database::Check>( data2 );
}
#endif


/************************************************************************
 * read/write HDF5 (AMP::DatabaseBox)                                   *
 ***********************************************************************/
#ifdef AMP_USE_HDF5
template<>
void AMP::writeHDF5Array<AMP::DatabaseBox>( hid_t fid,
                                            const std::string_view &name,
                                            const AMP::Array<DatabaseBox> &x )
{
    AMP::Array<int> y( cat( AMP::ArraySize( 11 ), x.size() ) );
    y.fill( 0 );
    for ( size_t i = 0; i < x.length(); i++ ) {
        int ndim  = x( i ).dim();
        y( 0, i ) = ndim;
        for ( int j = 0; j < ndim; j++ )
            y( 1 + j, i ) = x( i ).lower( j );
        for ( int j = 0; j < ndim; j++ )
            y( 6 + j, i ) = x( i ).upper( j );
    }
    writeHDF5( fid, name, y );
}
template<>
void AMP::readHDF5Array<AMP::DatabaseBox>( hid_t fid,
                                           const std::string_view &name,
                                           AMP::Array<DatabaseBox> &x )
{
    AMP::Array<int> y;
    readHDF5( fid, name, y );
    x.resize( pop( y.size() ) );
    for ( size_t i = 0; i < x.length(); i++ )
        x( i ) = DatabaseBox( y( 0, i ), &y( 1, i ), &y( 6, i ) );
}
template<>
void AMP::writeHDF5Scalar<AMP::DatabaseBox>( hid_t fid,
                                             const std::string_view &name,
                                             const DatabaseBox &x )
{
    auto y = AMP::Array<AMP::DatabaseBox>::staticView( { 1 }, const_cast<DatabaseBox *>( &x ) );
    writeHDF5Array<AMP::DatabaseBox>( fid, name, y );
}
template<>
void AMP::readHDF5Scalar<AMP::DatabaseBox>( hid_t fid,
                                            const std::string_view &name,
                                            DatabaseBox &x )
{
    AMP::Array<AMP::DatabaseBox> y;
    readHDF5Array<AMP::DatabaseBox>( fid, name, y );
    AMP_ASSERT( y.size() == AMP::ArraySize( 1 ) );
    x = y( 0 );
}
template void
AMP::readHDF5<AMP::DatabaseBox>( hid_t, const std::string_view &, AMP::DatabaseBox & );
template void
AMP::writeHDF5<AMP::DatabaseBox>( hid_t, const std::string_view &, const AMP::DatabaseBox & );
template void AMP::readHDF5<AMP::Array<AMP::DatabaseBox>>( hid_t,
                                                           const std::string_view &,
                                                           AMP::Array<AMP::DatabaseBox> & );
template void AMP::writeHDF5<AMP::Array<AMP::DatabaseBox>>( hid_t,
                                                            const std::string_view &,
                                                            const AMP::Array<AMP::DatabaseBox> & );
#else
INSTANTIATE_HDF5( AMP::DatabaseBox );
#endif


/************************************************************************
 * Explicit instantiations                                              *
 ***********************************************************************/
INSTANTIATE_HDF5( AMP::Database );
INSTANTIATE_HDF5( AMP::Database::Check );
