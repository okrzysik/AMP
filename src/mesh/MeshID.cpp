#include "AMP/mesh/MeshID.h"
#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/mesh/MeshIterator.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/Utilities.hpp"
#include "AMP/utils/kdtree2.hpp"


/********************************************************
 * MeshID                                                *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::MeshID ) == sizeof( uint64_t ) );
template<>
hid_t AMP::IO::getHDF5datatype<AMP::Mesh::MeshID>()
{
    return getHDF5datatype<uint64_t>();
}
template<>
void AMP::IO::writeHDF5Array<AMP::Mesh::MeshID>( hid_t fid,
                                                 const std::string &name,
                                                 const AMP::Array<AMP::Mesh::MeshID> &data )
{
    AMP::Array<uint64_t> data2( data.size(), reinterpret_cast<const uint64_t *>( data.data() ) );
    writeHDF5Array<uint64_t>( fid, name, data2 );
}
template<>
void AMP::IO::readHDF5Array<AMP::Mesh::MeshID>( hid_t fid,
                                                const std::string &name,
                                                AMP::Array<AMP::Mesh::MeshID> &data )
{
    Array<uint64_t> data2;
    readHDF5Array<uint64_t>( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = AMP::Mesh::MeshID( data2( i ) );
}
template<>
void AMP::IO::writeHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                                  const std::string &name,
                                                  const AMP::Mesh::MeshID &data )
{
    writeHDF5Scalar<uint64_t>( fid, name, data.getData() );
}
template<>
void AMP::IO::readHDF5Scalar<AMP::Mesh::MeshID>( hid_t fid,
                                                 const std::string &name,
                                                 AMP::Mesh::MeshID &data )
{
    uint64_t data2;
    readHDF5Scalar<uint64_t>( fid, name, data2 );
    data = AMP::Mesh::MeshID( data2 );
}
#endif


/********************************************************
 * MeshElementID                                         *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::MeshElementID ) == 2 * sizeof( uint64_t ) );
template<>
hid_t AMP::IO::getHDF5datatype<AMP::Mesh::MeshElementID>()
{
    return getHDF5datatype<uint64_t>();
}
template<>
void AMP::IO::writeHDF5Array<AMP::Mesh::MeshElementID>(
    hid_t fid, const std::string &name, const AMP::Array<AMP::Mesh::MeshElementID> &x )
{
    auto size2 = cat( ArraySize( 2 ), x.size() );
    auto ptr   = const_cast<uint64_t *>( reinterpret_cast<const uint64_t *>( x.data() ) );
    auto y     = AMP::Array<uint64_t>::staticView( size2, ptr );
    writeHDF5Array( fid, name, y );
}
template<>
void AMP::IO::readHDF5Array<AMP::Mesh::MeshElementID>( hid_t fid,
                                                       const std::string &name,
                                                       AMP::Array<AMP::Mesh::MeshElementID> &x )
{
    Array<uint64_t> y;
    readHDF5Array( fid, name, y );
    x.resize( pop( y.size() ) );
    for ( size_t i = 0; i < x.length(); i++ )
        x( i ) = AMP::Mesh::MeshElementID( AMP::Mesh::MeshID( y( 0, i ) ),
                                           AMP::Mesh::ElementID( y( 1, i ) ) );
}
template<>
void AMP::IO::writeHDF5Scalar<AMP::Mesh::MeshElementID>( hid_t fid,
                                                         const std::string &name,
                                                         const AMP::Mesh::MeshElementID &data )
{
    AMP::Array<AMP::Mesh::MeshElementID> x( { 1 }, &data );
    writeHDF5Array( fid, name, x );
}
template<>
void AMP::IO::readHDF5Scalar<AMP::Mesh::MeshElementID>( hid_t fid,
                                                        const std::string &name,
                                                        AMP::Mesh::MeshElementID &data )
{
    AMP::Array<AMP::Mesh::MeshElementID> x;
    readHDF5Array( fid, name, x );
    AMP_ASSERT( x.size() == 1u );
    data = x( 0 );
}
#endif
INSTANTIATE_HDF5( AMP::Mesh::MeshElementID );


/********************************************************
 * GeomType                                              *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::GeomType ) == sizeof( uint8_t ) );
template<>
hid_t AMP::IO::getHDF5datatype<AMP::Mesh::GeomType>()
{
    return getHDF5datatype<uint8_t>();
}
template<>
void AMP::IO::writeHDF5Array<AMP::Mesh::GeomType>( hid_t fid,
                                                   const std::string &name,
                                                   const AMP::Array<AMP::Mesh::GeomType> &data )
{
    Array<uint8_t> data2( data.size(), reinterpret_cast<const uint8_t *>( data.data() ) );
    writeHDF5Array<uint8_t>( fid, name, data2 );
}
template<>
void AMP::IO::readHDF5Array<AMP::Mesh::GeomType>( hid_t fid,
                                                  const std::string &name,
                                                  AMP::Array<AMP::Mesh::GeomType> &data )
{
    Array<uint8_t> data2;
    readHDF5Array<uint8_t>( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = static_cast<AMP::Mesh::GeomType>( data2( i ) );
}
template<>
void AMP::IO::writeHDF5Scalar<AMP::Mesh::GeomType>( hid_t fid,
                                                    const std::string &name,
                                                    const AMP::Mesh::GeomType &data )
{
    writeHDF5Scalar<uint8_t>( fid, name, static_cast<uint8_t>( data ) );
}
template<>
void AMP::IO::readHDF5Scalar<AMP::Mesh::GeomType>( hid_t fid,
                                                   const std::string &name,
                                                   AMP::Mesh::GeomType &data )
{
    uint8_t data2;
    readHDF5Scalar<uint8_t>( fid, name, data2 );
    data = static_cast<AMP::Mesh::GeomType>( data2 );
}
#endif


/********************************************************
 * AMP::Mesh::MeshIterator::Type                         *
 ********************************************************/
#ifdef AMP_USE_HDF5
static_assert( sizeof( AMP::Mesh::MeshIterator::Type ) == sizeof( uint8_t ) );
template<>
hid_t AMP::IO::getHDF5datatype<AMP::Mesh::MeshIterator::Type>()
{
    return getHDF5datatype<uint8_t>();
}
template<>
void AMP::IO::writeHDF5Array<AMP::Mesh::MeshIterator::Type>(
    hid_t fid, const std::string &name, const AMP::Array<AMP::Mesh::MeshIterator::Type> &data )
{
    Array<uint8_t> data2( data.size(), reinterpret_cast<const uint8_t *>( data.data() ) );
    writeHDF5Array<uint8_t>( fid, name, data2 );
}
template<>
void AMP::IO::readHDF5Array<AMP::Mesh::MeshIterator::Type>(
    hid_t fid, const std::string &name, AMP::Array<AMP::Mesh::MeshIterator::Type> &data )
{
    Array<uint8_t> data2;
    readHDF5Array<uint8_t>( fid, name, data2 );
    data.resize( data2.size() );
    for ( size_t i = 0; i < data.length(); i++ )
        data( i ) = static_cast<AMP::Mesh::MeshIterator::Type>( data2( i ) );
}
template<>
void AMP::IO::writeHDF5Scalar<AMP::Mesh::MeshIterator::Type>(
    hid_t fid, const std::string &name, const AMP::Mesh::MeshIterator::Type &data )
{
    writeHDF5Scalar<uint8_t>( fid, name, static_cast<uint8_t>( data ) );
}
template<>
void AMP::IO::readHDF5Scalar<AMP::Mesh::MeshIterator::Type>( hid_t fid,
                                                             const std::string &name,
                                                             AMP::Mesh::MeshIterator::Type &data )
{
    uint8_t data2;
    readHDF5Scalar<uint8_t>( fid, name, data2 );
    data = static_cast<AMP::Mesh::MeshIterator::Type>( data2 );
}
#endif
INSTANTIATE_HDF5( AMP::Mesh::MeshIterator::Type );


/********************************************************
 * Explicit instantiations                               *
 ********************************************************/
#define INSTANTIATE_SORT( TYPE )                                                   \
    template void AMP::Utilities::quicksort<TYPE>( size_t, TYPE * );               \
    template void AMP::Utilities::quicksort<TYPE, TYPE>( size_t, TYPE *, TYPE * ); \
    template void AMP::Utilities::unique<TYPE>( std::vector<TYPE> & );             \
    template void AMP::Utilities::unique<TYPE>(                                    \
        std::vector<TYPE> &, std::vector<size_t> &, std::vector<size_t> & );       \
    template size_t AMP::Utilities::findfirst<TYPE>( size_t, const TYPE *, const TYPE & )
using ElementArray1 = std::array<AMP::Mesh::ElementID, 1>;
using ElementArray2 = std::array<AMP::Mesh::ElementID, 2>;
using ElementArray3 = std::array<AMP::Mesh::ElementID, 3>;
using ElementArray4 = std::array<AMP::Mesh::ElementID, 4>;
INSTANTIATE_SORT( AMP::Mesh::GeomType );
INSTANTIATE_SORT( AMP::Mesh::MeshID );
INSTANTIATE_SORT( AMP::Mesh::ElementID );
INSTANTIATE_SORT( AMP::Mesh::MeshElementID );
INSTANTIATE_SORT( ElementArray1 );
INSTANTIATE_SORT( ElementArray2 );
INSTANTIATE_SORT( ElementArray3 );
INSTANTIATE_SORT( ElementArray4 );
INSTANTIATE_HDF5( AMP::Mesh::GeomType );
INSTANTIATE_HDF5( AMP::Mesh::MeshID );
instantiateArrayConstructors( AMP::Mesh::MeshID );
instantiateArrayConstructors( AMP::Mesh::GeomType );
template class AMP::kdtree2<1, AMP::Mesh::MeshElementID>;
template class AMP::kdtree2<2, AMP::Mesh::MeshElementID>;
template class AMP::kdtree2<3, AMP::Mesh::MeshElementID>;
template void AMP::Utilities::quicksort<int, AMP::Mesh::MeshElementID>(
    size_t, int *, AMP::Mesh::MeshElementID * );
