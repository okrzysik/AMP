#include "AMP/mesh/MeshID.h"
#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/utils/Array.hpp"


/********************************************************
 * MeshID                                                *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::MeshID ) == sizeof( uint64_t ) );
template<>
hid_t AMP::getHDF5datatype<AMP::Mesh::MeshID>()
{
    return getHDF5datatype<uint64_t>();
}
template<>
void AMP::writeHDF5Array<AMP::Mesh::MeshID>( hid_t fid,
                                             const std::string_view &name,
                                             const AMP::Array<AMP::Mesh::MeshID> &data )
{
    AMP::Array<uint64_t> data2( data.size(), reinterpret_cast<const uint64_t *>( data.data() ) );
    writeHDF5Array<uint64_t>( fid, name, data2 );
}
template<>
void AMP::readHDF5Array<AMP::Mesh::MeshID>( hid_t fid,
                                            const std::string_view &name,
                                            AMP::Array<AMP::Mesh::MeshID> &data )
{
    AMP::Array<uint64_t> data2;
    AMP::readHDF5Array<uint64_t>( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = AMP::Mesh::MeshID( data2( i ) );
}
template<>
void AMP::writeHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                              const std::string_view &name,
                                              const AMP::Mesh::MeshID &data )
{
    writeHDF5Scalar<uint64_t>( fid, name, data.getData() );
}
template<>
void AMP::readHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                             const std::string_view &name,
                                             AMP::Mesh::MeshID &data )
{
    uint64_t data2;
    readHDF5Scalar<uint64_t>( fid, name, data2 );
    data = AMP::Mesh::MeshID( data2 );
}
#endif


/********************************************************
 * GeomType                                              *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::GeomType ) == sizeof( uint8_t ) );
template<>
hid_t AMP::getHDF5datatype<AMP::Mesh::GeomType>()
{
    return getHDF5datatype<uint16_t>();
}
template<>
void AMP::writeHDF5Array<AMP::Mesh::GeomType>( hid_t fid,
                                               const std::string_view &name,
                                               const AMP::Array<AMP::Mesh::GeomType> &data )
{
    AMP::Array<uint16_t> data2( data.size(), reinterpret_cast<const uint16_t *>( data.data() ) );
    writeHDF5Array<uint16_t>( fid, name, data2 );
}
template<>
void AMP::readHDF5Array<AMP::Mesh::GeomType>( hid_t fid,
                                              const std::string_view &name,
                                              AMP::Array<AMP::Mesh::GeomType> &data )
{
    AMP::Array<uint16_t> data2;
    AMP::readHDF5Array<uint16_t>( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = static_cast<AMP::Mesh::GeomType>( data2( i ) );
}
template<>
void AMP::writeHDF5Scalar<AMP::Mesh::GeomType>( hid_t fid,
                                                const std::string_view &name,
                                                const AMP::Mesh::GeomType &data )
{
    writeHDF5Scalar<uint16_t>( fid, name, static_cast<uint16_t>( data ) );
}
template<>
void AMP::readHDF5Scalar<AMP::Mesh::GeomType>( hid_t fid,
                                               const std::string_view &name,
                                               AMP::Mesh::GeomType &data )
{
    uint16_t data2;
    readHDF5Scalar<uint16_t>( fid, name, data2 );
    data = static_cast<AMP::Mesh::GeomType>( data2 );
}
#endif


/********************************************************
 * Explicit instantiations                               *
 ********************************************************/
INSTANTIATE_HDF5( AMP::Mesh::GeomType );
INSTANTIATE_HDF5( AMP::Mesh::MeshID );
instantiateArrayConstructors( AMP::Mesh::MeshID );
instantiateArrayConstructors( AMP::Mesh::GeomType );
