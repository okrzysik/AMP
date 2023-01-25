#include "AMP/mesh/MeshID.h"
#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/utils/Array.hpp"


/********************************************************
 * HDF5 operators                                        *
 ********************************************************/
#ifdef AMP_USE_HDF5
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
    writeHDF5( fid, name, data2 );
}
template<>
void AMP::readHDF5Array<AMP::Mesh::MeshID>( hid_t fid,
                                            const std::string_view &name,
                                            AMP::Array<AMP::Mesh::MeshID> &data )
{
    AMP::Array<uint64_t> data2;
    AMP::readHDF5( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = AMP::Mesh::MeshID( data2( i ) );
}
template<>
void AMP::writeHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                              const std::string_view &name,
                                              const AMP::Mesh::MeshID &data )
{
    writeHDF5( fid, name, data.getData() );
}
template<>
void AMP::readHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                             const std::string_view &name,
                                             AMP::Mesh::MeshID &data )
{
    uint64_t data2;
    readHDF5( fid, name, data2 );
    data = AMP::Mesh::MeshID( data2 );
}
INSTANTIATE_HDF5( AMP::Mesh::MeshID );
#endif


/********************************************************
 * HDF5 operators                                        *
 ********************************************************/
instantiateArrayConstructors( AMP::Mesh::MeshID );
