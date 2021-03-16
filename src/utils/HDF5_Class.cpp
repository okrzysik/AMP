#include "AMP/utils/HDF5_Class.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/HDF5_IO.h"
#include "AMP/utils/HDF5_IO.hpp"
#include "AMP/utils/PIO.h"
#include "AMP/utils/Utilities.h"

#include <set>
#include <sstream>
#include <utility>
#include <vector>


#ifdef USE_HDF5 // USE HDF5


namespace AMP {


static inline std::string HDF5_getMemberName( hid_t id, unsigned idx )
{
    char *cname = H5Tget_member_name( id, idx );
    std::string name( cname );
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 8
    free( cname );
#else
    H5free_memory( cname );
#endif
    return name;
}


/******************************************************************
 * Classes to store HDF5 data                                      *
 ******************************************************************/
static int find( const std::vector<std::string> &vec, const std::string_view &x )
{
    for ( size_t i = 0; i < vec.size(); i++ ) {
        if ( vec[i] == x )
            return i;
    }
    return -1;
}
class HDF5_null final : public HDF5data
{
public:
    HDF5_null( hid_t fid, const std::string_view &name, const std::string_view &type )
        : HDF5data( fid, name ), d_type( std::move( type ) )
    {
    }
    ~HDF5_null() override = default;
    std::string type() const override { return "HDF5_null"; }
    size_t size() const override { return 0; }
    std::shared_ptr<HDF5data> getData( size_t, const std::string_view & ) override
    {
        return nullptr;
    }
    std::vector<std::string> getNames() const override { return std::vector<std::string>(); }
    void print( int, const std::string_view &prefix ) const override
    {
        printf( "%s%s - Unfinished %s\n", prefix.data(), d_name.data(), d_type.data() );
    }

private:
    std::string d_type;
};
class HDF5_group final : public HDF5data
{
public:
    HDF5_group( hid_t fid, const std::string_view &name, int type );
    ~HDF5_group() override = default;
    std::string type() const override { return "HDF5_group"; }
    size_t size() const override { return d_data.length() / d_data.size( 0 ); }
    std::shared_ptr<HDF5data> getData( size_t i, const std::string_view &name ) override
    {
        int j = find( d_names, name );
        if ( j == -1 )
            return nullptr;
        return d_data( j, i );
    }
    std::vector<std::string> getNames() const override { return d_names; }
    void print( int level, const std::string_view &prefix = "" ) const override
    {
        if ( d_data.empty() )
            return;
        size_t N = d_data.length() / d_data.size( 0 );
        if ( N == 1 && level >= 2 ) {
            printf( "%s%s\n", prefix.data(), d_name.data() );
            for ( size_t i = 0; i < d_data.size( 0 ); i++ ) {
                if ( d_data( i ) )
                    d_data( i )->print( level, std::string( prefix ) + "  " );
            }
        } else {
            printf( "%s%s (%i", prefix.data(), d_name.data(), (int) d_data.size( 1 ) );
            for ( int d = 2; d < d_data.ndim(); d++ )
                printf( ",%i", (int) d_data.size( d ) );
            printf( ")\n" );
            size_t N_max = N;
            if ( level < 3 )
                N_max = std::min<size_t>( N_max, 5 );
            for ( size_t j = 0; j < N_max; j++ ) {
                printf( "%s  [%i]\n", prefix.data(), static_cast<int>( j + 1 ) );
                for ( size_t i = 0; i < d_data.size( 0 ); i++ ) {
                    if ( d_data( i, j ) && level >= 2 )
                        d_data( i, j )->print( level, std::string( prefix ) + "    " );
                }
            }
            if ( N > N_max )
                printf( "%s  ...\n", prefix.data() );
        }
    }

private:
    // Store pointers to the data
    // Note: the first dimension is the number of variables
    //    All other dimensions are the dimensions of this
    std::vector<std::string> d_names;
    AMP::Array<std::shared_ptr<HDF5data>> d_data;
};
template<class TYPE>
class HDF5_primitive final : public HDF5data
{
public:
    HDF5_primitive( hid_t fid, const std::string_view &name );
    HDF5_primitive( const std::string_view &name, const AMP::Array<TYPE> &data );
    ~HDF5_primitive() override = default;
    std::string type() const override
    {
        return AMP::Utilities::stringf( "HDF5_primitive<%f>", typeid( TYPE ).name() );
    }
    size_t size() const override { return 1; }
    std::shared_ptr<HDF5data> getData( size_t i, const std::string_view &name ) override
    {
        if ( i == 0 && name == d_name )
            return shared_from_this();
        return nullptr;
    }
    std::vector<std::string> getNames() const override
    {
        return std::vector<std::string>( 1, d_name );
    }
    const AMP::Array<TYPE> &getData() const { return d_data; }
    void print( int level, const std::string_view &prefix = "" ) const override
    {
        printf( "%s%s (%s)", prefix.data(), d_name.data(), typeid( TYPE ).name() );
        if ( d_data.empty() ) {
            printf( " []\n" );
        } else if ( d_data.length() == 1 ) {
            AMP::pout << ": " << d_data( 0 ) << std::endl;
        } else if ( d_data.length() <= 4 || level >= 4 ) {
            AMP::pout << " [";
            if constexpr ( std::is_same<TYPE, unsigned char>::value ) {
                AMP::pout << " " << static_cast<int>( d_data( 0 ) );
                for ( size_t i = 1; i < d_data.length(); i++ )
                    AMP::pout << ", " << static_cast<int>( d_data( i ) );
            } else {
                AMP::pout << " " << d_data( 0 );
                for ( size_t i = 1; i < d_data.length(); i++ )
                    AMP::pout << ", " << d_data( i );
            }
            AMP::pout << " ]" << std::endl;
        } else {
            printf( " (%i", (int) d_data.size( 0 ) );
            for ( int d = 1; d < d_data.ndim(); d++ )
                printf( ",%i", (int) d_data.size( d ) );
            printf( ")\n" );
        }
    }

private:
    AMP::Array<TYPE> d_data;
};
template<>
void HDF5data::getData<std::string>( AMP::Array<std::string> &data ) const
{
    auto tmp1 = dynamic_cast<const HDF5_primitive<std::string> *>( this );
    auto tmp2 = dynamic_cast<const HDF5_primitive<char> *>( this );
    if ( tmp1 ) {
        data = tmp1->getData();
    } else if ( tmp2 ) {
        auto tmp = tmp2->getData();
        data.resize( 1 );
        data( 0 ) = std::string( tmp.data(), tmp.length() );
    } else {
        AMP_ERROR( "Internal error" );
    }
}
template<class TYPE>
void HDF5data::getData( AMP::Array<TYPE> &data ) const
{
    if ( dynamic_cast<const HDF5_primitive<char> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<char> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<unsigned char> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<unsigned char> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<int> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<int> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<unsigned int> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<unsigned int> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<long int> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<long int> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<unsigned long int> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<unsigned long int> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<float> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<float> *>( this )->getData() );
    } else if ( dynamic_cast<const HDF5_primitive<double> *>( this ) != nullptr ) {
        data.copy( dynamic_cast<const HDF5_primitive<double> *>( this )->getData() );
    } else {
        AMP_ERROR( "Unable to get data: " + type() );
    }
}
template void HDF5data::getData<int>( AMP::Array<int> & ) const;
template void HDF5data::getData<char>( AMP::Array<char> & ) const;
template void HDF5data::getData<long>( AMP::Array<long> & ) const;
template void HDF5data::getData<float>( AMP::Array<float> & ) const;
template void HDF5data::getData<double>( AMP::Array<double> & ) const;
template void HDF5data::getData<unsigned int>( AMP::Array<unsigned int> & ) const;
template void HDF5data::getData<unsigned char>( AMP::Array<unsigned char> & ) const;
template void HDF5data::getData<unsigned long>( AMP::Array<unsigned long> & ) const;


/************************************************************************
 * Read database entry                                                   *
 ************************************************************************/
template<class TYPE>
HDF5_primitive<TYPE>::HDF5_primitive( hid_t fid, const std::string_view &name )
    : HDF5data( fid, name )
{
    readHDF5( fid, name, d_data );
}
template<class TYPE>
HDF5_primitive<TYPE>::HDF5_primitive( const std::string_view &name, const AMP::Array<TYPE> &data )
    : HDF5data( 0, name ), d_data( std::move( data ) )
{
}
static std::unique_ptr<HDF5data> readPrimitive( hid_t fid, const std::string_view &name )
{
    hid_t id  = H5Dopen( fid, name.data(), H5P_DEFAULT );
    hid_t tid = H5Dget_type( id );
    std::unique_ptr<HDF5data> data;
    if ( H5Tequal( tid, getHDF5datatype<char>() ) ) {
        data.reset( new HDF5_primitive<char>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<unsigned char>() ) ) {
        data.reset( new HDF5_primitive<unsigned char>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<int>() ) ) {
        data.reset( new HDF5_primitive<int>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<unsigned int>() ) ) {
        data.reset( new HDF5_primitive<unsigned int>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<long int>() ) ) {
        data.reset( new HDF5_primitive<long int>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<unsigned long int>() ) ) {
        data.reset( new HDF5_primitive<unsigned long int>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<float>() ) ) {
        data.reset( new HDF5_primitive<float>( fid, name ) );
    } else if ( H5Tequal( tid, getHDF5datatype<double>() ) ) {
        data.reset( new HDF5_primitive<double>( fid, name ) );
    } else {
        AMP_ERROR( "Unknown data" );
    }
    return data;
}
static std::unique_ptr<HDF5data> readDatabase( hid_t fid, const std::string_view &name )
{
    hid_t id            = H5Dopen( fid, name.data(), H5P_DEFAULT );
    hid_t tid           = H5Dget_type( id );
    H5T_class_t classid = H5Tget_class( tid );
    std::unique_ptr<HDF5data> data;
    if ( classid == H5T_INTEGER || classid == H5T_FLOAT ) {
        data = readPrimitive( fid, name );
    } else if ( classid == H5T_STRING ) {
        data.reset( new HDF5_primitive<std::string>( fid, name ) );
    } else if ( classid == H5T_BITFIELD ) {
        data.reset( new HDF5_null( fid, name, "H5T_BITFIELD" ) );
    } else if ( classid == H5T_OPAQUE ) {
        data.reset( new HDF5_null( fid, name, "H5T_OPAQUE" ) );
    } else if ( classid == H5T_REFERENCE ) {
        data.reset( new HDF5_null( fid, name, "H5T_REFERENCE" ) );
    } else if ( classid == H5T_ENUM ) {
        data.reset( new HDF5_null( fid, name, "H5T_ENUM" ) );
    } else if ( classid == H5T_VLEN ) {
        data.reset( new HDF5_null( fid, name, "H5T_VLEN" ) );
    } else if ( classid == H5T_ARRAY ) {
        data.reset( new HDF5_null( fid, name, "H5T_ARRAY" ) );
    } else if ( classid == H5T_COMPOUND ) {
        data.reset( new HDF5_group( fid, name, 2 ) );
    } else {
        AMP_ERROR( "Unknown data" );
    }
    return data;
}


/************************************************************************
 * Read group / compound data                                            *
 ************************************************************************/
template<class TYPE>
static void getGroupData( const std::string_view &name,
                          size_t size,
                          const char *data,
                          size_t N,
                          size_t offset,
                          size_t varnum,
                          AMP::Array<std::shared_ptr<HDF5data>> &out )
{
    for ( size_t i = 0; i < N; i++ ) {
        AMP::Array<TYPE> var( 1 );
        memcpy( var.data(), &data[i * size + offset], sizeof( TYPE ) );
        out( varnum, i ).reset( new HDF5_primitive<TYPE>( name, var ) );
    }
}
size_t getIndex( const std::string_view &name )
{
    auto i = name.find_first_of( "0123456789" );
    if ( i == std::string::npos )
        return 0;
    auto substr = name.substr( i );
    auto index  = std::stoi( substr.data() );
    return index;
}
HDF5_group::HDF5_group( hid_t fid, const std::string_view &name, int type ) : HDF5data( fid, name )
{
    if ( type == 1 ) {
        // Read group
        hid_t gid = H5Gopen2( fid, name.data(), H5P_DEFAULT );
        H5G_info_t group_info;
        H5Gget_info( gid, &group_info );
        std::vector<std::shared_ptr<HDF5data>> data;
        data.reserve( group_info.nlinks );
        for ( size_t i = 0; i < group_info.nlinks; i++ ) {
            char name2[512];
            H5Gget_objname_by_idx( gid, i, name2, sizeof( name2 ) );
            if ( std::string( name2 ) == ".." )
                continue;
            auto data2 = readHDF5( gid, name2 );
            if ( data2 )
                data.push_back( std::move( data2 ) );
            d_names.emplace_back( name2 );
            AMP_ASSERT( data.back() );
        }
        d_data = data;
        // Check if variables represent an array
        bool check = true;
        int N      = -1;
        int imin   = 100000;
        int imax   = -1;
        std::vector<bool> test( d_names.size(), false );
        for ( size_t i = 0; i < d_names.size(); i++ ) {
            bool containsLetter =
                d_names[i].find_first_of(
                    "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" ) != std::string::npos;
            if ( d_names[i] == "N" ||
                 d_names[i].compare( 0, 2 + d_name.size(), "N_" + d_name ) == 0 ) {
                AMP::Array<int> data;
                d_data( i )->getData( data );
                AMP_ASSERT( data.length() == 1u );
                N       = std::max<int>( N, data( 0 ) );
                test[i] = false;
            } else if ( d_names[i].compare( 0, d_name.size(), d_name ) == 0 ) {
                size_t index = getIndex( d_names[i] );
                imin         = std::min<int>( imin, index );
                imax         = std::max<int>( imax, index );
                test[i]      = true;
            } else if ( !containsLetter ) {
                size_t index = getIndex( d_names[i] );
                imin         = std::min<int>( imin, index );
                imax         = std::max<int>( imax, index );
                test[i]      = true;
            } else {
                check   = false;
                test[i] = false;
            }
        }
        int offset = imin == 0 ? 0 : -1;
        if ( N == -1 )
            N = imax + offset + 1;
        check = check && ( N + 1 ) >= (int) d_data.length() && imax + offset < N;
        if ( check ) {
            // Collapse the variables into an array
            AMP_ASSERT( d_data.length() == d_names.size() );
            auto old_vars = d_names;
            auto old_data = d_data;
            // Get a list of the new variabes
            std::set<std::string> set;
            for ( const auto &tmp : old_data ) {
                for ( const auto &tmp2 : tmp->getNames() ) {
                    set.insert( tmp2 );
                }
            }
            d_names = std::vector<std::string>( set.begin(), set.end() );
            // Create the new structure
            d_data = AMP::Array<std::shared_ptr<HDF5data>>( d_names.size(), N );
            for ( size_t i = 0; i < old_vars.size(); i++ ) {
                if ( test[i] ) {
                    size_t index = getIndex( old_vars[i] ) + offset;
                    AMP_ASSERT( old_data( i )->size() == 1u );
                    for ( size_t j = 0; j < d_names.size(); j++ )
                        d_data( j, index ) = old_data( i )->getData( 0, d_names[j] );
                }
            }
        }
        H5Gclose( gid );
    } else if ( type == 2 ) {
        // Read compound array as a group
        hid_t id        = H5Dopen2( fid, name.data(), H5P_DEFAULT );
        hid_t tid       = H5Dget_type( id );
        hid_t dataspace = H5Dget_space( id );
        hsize_t dims0[10];
        int ndim  = H5Sget_simple_extent_dims( dataspace, dims0, nullptr );
        auto dims = convertSize( ndim, dims0 );
        size_t N  = 1;
        for ( auto dim : dims )
            N *= dim;
        size_t sizeofobj = H5Tget_size( tid );
        auto data        = new char[sizeofobj * N];
        memset( data, 0, sizeofobj * N );
        H5Dread( id, tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
        int N_members = H5Tget_nmembers( tid );
        std::vector<size_t> dims2( 1, N_members );
        for ( auto dim : dims )
            dims2.emplace_back( dim );
        d_data.resize( dims2 );
        for ( int i = 0; i < N_members; i++ ) {
            hid_t tid2    = H5Tget_member_type( tid, i );
            auto name2    = HDF5_getMemberName( tid, i );
            size_t offset = H5Tget_member_offset( tid, i );
            d_names.emplace_back( name2 );
            if ( H5Tequal( tid2, getHDF5datatype<char>() ) ) {
                getGroupData<char>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<unsigned char>() ) ) {
                getGroupData<unsigned char>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<int>() ) ) {
                getGroupData<int>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<unsigned int>() ) ) {
                getGroupData<unsigned int>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<long int>() ) ) {
                getGroupData<long int>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<unsigned long int>() ) ) {
                getGroupData<unsigned long int>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<float>() ) ) {
                getGroupData<float>( name2, sizeofobj, data, N, offset, i, d_data );
            } else if ( H5Tequal( tid2, getHDF5datatype<double>() ) ) {
                getGroupData<double>( name2, sizeofobj, data, N, offset, i, d_data );
            } else {
                AMP_ERROR( "Unknown data" );
            }
        }
        delete[] data;
        H5Dclose( id );
        H5Tclose( tid );
        H5Sclose( dataspace );
    }
}


/************************************************************************
 * Read arbitrary data                                                   *
 ************************************************************************/
std::unique_ptr<HDF5data> readHDF5( hid_t fid, const std::string_view &name )
{
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR <= 8
    H5O_info_t object_info;
    H5Oget_info_by_name( fid, name.data(), &object_info, H5P_DEFAULT );
#else
    H5O_info1_t object_info;
    H5Oget_info_by_name1( fid, name.data(), &object_info, H5P_DEFAULT );
#endif
    std::unique_ptr<HDF5data> data;
    if ( object_info.type == H5O_TYPE_GROUP ) {
        data.reset( new HDF5_group( fid, name, 1 ) );
    } else if ( object_info.type == H5O_TYPE_DATASET ) {
        data = readDatabase( fid, name );
    } else if ( object_info.type == H5O_TYPE_NAMED_DATATYPE ) {
        AMP_WARNING( "Read named datatype not finished: " + std::string( name ) );
    } else if ( object_info.type == H5O_TYPE_UNKNOWN ) {
        AMP_ERROR( "Read unknown type not supported" + std::string( name ) );
    } else {
        AMP_ERROR( "Read Error" );
    }
    return data;
}


} // namespace AMP

#endif
