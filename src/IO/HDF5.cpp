#include "AMP/IO/HDF5.h"
#include "AMP/IO/HDF5.hpp"
#include "AMP/utils/Array.h"
#include "AMP/utils/Utilities.h"

#include <complex>
#include <sstream>
#include <string>
#include <vector>


namespace AMP {


#ifdef AMP_USE_HDF5 // USE HDF5


/******************************************************************
 * Open/close HDF5 files                                           *
 ******************************************************************/
hid_t openHDF5( const std::string_view &filename, const char *mode, Compression compress )
{
    // Set cache size to 3MBs and instruct the cache to discard the fully read chunk
    auto pid = H5P_DEFAULT;
    /*auto pid = H5Pcreate( H5P_FILE_ACCESS );
    int nelemts;
    size_t nslots, nbytes;
    double w0;
    H5Pget_cache(pid,& nelemts,& nslots,& nbytes,& w0);
    H5Pset_cache(pid, nelemts, 1999, 3*1024*1024, 1.0); */
    // Open the file
    hid_t fid = 0;
    if ( strcmp( mode, "r" ) == 0 ) {
        fid = H5Fopen( filename.data(), H5F_ACC_RDONLY, pid );
    } else if ( strcmp( mode, "w" ) == 0 ) {
        fid = H5Fcreate( filename.data(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    } else if ( strcmp( mode, "rw" ) == 0 ) {
        fid = H5Fopen( filename.data(), H5F_ACC_RDWR, H5P_DEFAULT );
    } else {
        AMP_ERROR( "Invalid mode for opening HDF5 file" );
    }
    if ( strcmp( mode, "w" ) == 0 ) {
        int defaultCompression = 0;
        if ( compress == Compression::None ) {
            // No compression
        } else if ( compress == Compression::GZIP ) {
            // Use gzip (if available)
            if ( H5Zfilter_avail( H5Z_FILTER_DEFLATE ) )
                defaultCompression = 1;
            else
                AMP::perr << "HDF5 gzip filter is not available" << std::endl;
        } else if ( compress == Compression::SZIP ) {
            // Use szip (if available)
            if ( H5Zfilter_avail( H5Z_FILTER_SZIP ) )
                defaultCompression = 2;
            else
                AMP::perr << "HDF5 szip filter is not available" << std::endl;
        } else {
            AMP_ERROR( "Internal error" );
        }
        writeHDF5<int>( fid, "DefaultCompression", defaultCompression );
    }
    // H5Pclose( pid );
    return fid;
}
void closeHDF5( hid_t fid )
{
    // Try to close any remaining objects (needed to ensure we can reopen the data if desired)
    hid_t file[1000], set[1000], group[1000], type[1000], attr[1000];
    size_t N_file  = H5Fget_obj_ids( fid, H5F_OBJ_FILE, 1000, file );
    size_t N_set   = H5Fget_obj_ids( fid, H5F_OBJ_DATASET, 1000, set );
    size_t N_group = H5Fget_obj_ids( fid, H5F_OBJ_GROUP, 1000, group );
    size_t N_type  = H5Fget_obj_ids( fid, H5F_OBJ_DATATYPE, 1000, type );
    size_t N_attr  = H5Fget_obj_ids( fid, H5F_OBJ_ATTR, 1000, attr );
    for ( size_t i = 0; i < N_file; i++ ) {
        if ( file[i] != fid )
            H5Fclose( file[i] );
    }
    for ( size_t i = 0; i < N_set; i++ )
        H5Dclose( set[i] );
    for ( size_t i = 0; i < N_group; i++ )
        H5Gclose( group[i] );
    for ( size_t i = 0; i < N_type; i++ )
        H5Tclose( type[i] );
    for ( size_t i = 0; i < N_attr; i++ )
        H5Aclose( attr[i] );
    // Flush the data (needed to ensure we can reopen the data if desired)
    unsigned intent;
    H5Fget_intent( fid, &intent );
    if ( intent == H5F_ACC_RDWR || intent == H5F_ACC_TRUNC )
        H5Fflush( fid, H5F_SCOPE_GLOBAL );
    // Close the file
    H5Fclose( fid );
}


/************************************************************************
 * Check if we support compression                                       *
 ************************************************************************/
Compression defaultCompression( hid_t fid )
{
    hid_t root = H5Gopen2( fid, "/", H5P_DEFAULT );
    if ( !H5Dexists( root, "DefaultCompression" ) )
        return Compression::None;
    int tmp;
    readHDF5( root, "DefaultCompression", tmp );
    Compression compress = Compression::None;
    if ( tmp == 0 ) {
        compress = Compression::None;
    } else if ( tmp == 1 ) {
        compress = Compression::GZIP;
    } else if ( tmp == 2 ) {
        compress = Compression::SZIP;
    } else {
        AMP_ERROR( "Internal error" );
    }
    return compress;
}


/************************************************************************
 * Create a default chunk size                                           *
 ************************************************************************/
static size_t ceil2( size_t x )
{
    size_t y = 2;
    while ( y < x )
        y <<= 1;
    return y;
}
static size_t floor2( size_t x ) { return ceil2( x ) >> 1; }
hid_t createChunk( AMP::ArraySize dims, Compression compress, size_t objSize )
{
    if ( compress == Compression::None || dims.length() < 512 )
        return H5P_DEFAULT;
    if ( objSize == 0 )
        objSize = 512;
    // Limit chunk size to < 2 GB
    int d = dims.ndim() - 1;
    while ( dims.length() * objSize >= 0x80000000 ) {
        if ( dims[d] == 1 ) {
            d--;
        } else {
            dims.resize( d, floor2( dims[d] ) );
        }
    }
    hid_t plist = H5Pcreate( H5P_DATASET_CREATE );
    auto status = H5Pset_chunk( plist, dims.ndim(), arraySize( dims ).data() );
    AMP_ASSERT( status == 0 );
    if ( compress == Compression::GZIP ) {
        status = H5Pset_deflate( plist, 7 );
        AMP_ASSERT( status == 0 );
    } else if ( compress == Compression::SZIP ) {
        status = H5Pset_szip( plist, H5_SZIP_NN_OPTION_MASK, 16 );
        AMP_ASSERT( status == 0 );
    }
    return plist;
}


/************************************************************************
 * HDF5 helper routines                                                  *
 ************************************************************************/
bool H5Gexists( hid_t fid, const std::string_view &name )
{
    H5E_auto2_t func;
    void *client;
    H5Eget_auto2( H5E_DEFAULT, &func, &client );
    H5Eset_auto2( H5E_DEFAULT, nullptr, nullptr );
    int status = H5Gget_objinfo( fid, name.data(), 0, nullptr );
    H5Eset_auto2( H5E_DEFAULT, func, client );
    return status == 0;
}
bool H5Dexists( hid_t fid, const std::string_view &name )
{
    H5E_auto2_t func;
    void *client;
    H5Eget_auto2( H5E_DEFAULT, &func, &client );
    H5Eset_auto2( H5E_DEFAULT, nullptr, nullptr );
    hid_t dataset = H5Dopen2( fid, name.data(), H5P_DEFAULT );
    H5Eset_auto2( H5E_DEFAULT, func, client );
    bool exists = dataset > 0;
    // if ( exists )
    //    H5Dclose( dataset );
    return exists;
}
hid_t createGroup( hid_t fid, const std::string_view &name )
{
    return H5Gcreate2( fid, name.data(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
}
hid_t openGroup( hid_t fid, const std::string_view &name )
{
    if ( !H5Gexists( fid, name ) )
        AMP_ERROR( "Group " + std::string( name ) + " does not exist" );
    return H5Gopen2( fid, name.data(), H5P_DEFAULT );
}
void closeGroup( hid_t gid ) { H5Gclose( gid ); }


/************************************************************************
 * Read/write bytes to HDF5                                              *
 ************************************************************************/
void writeHDF5( hid_t fid, const std::string_view &name, size_t N_bytes, const void *data )
{
    hsize_t dim[1] = { N_bytes };
    hid_t plist    = H5P_DEFAULT;
    if ( N_bytes < 0x7500 ) {
        // Use compact storage (limited to < 30K)
        plist       = H5Pcreate( H5P_DATASET_CREATE );
        auto status = H5Pset_layout( plist, H5D_COMPACT );
        AMP_ASSERT( status == 0 );
    } else {
        // Use compression if available
        plist = createChunk( N_bytes, defaultCompression( fid ), 1 );
    }
    hid_t dataspace = H5Screate_simple( 1, dim, NULL );
    hid_t dataset   = H5Dcreate2(
        fid, name.data(), H5T_NATIVE_UCHAR, dataspace, H5P_DEFAULT, plist, H5P_DEFAULT );
    H5Dwrite( dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    H5Sclose( dataspace );
    if ( plist != H5P_DEFAULT )
        H5Pclose( plist );
}
void readHDF5( hid_t fid, const std::string_view &name, size_t N_bytes, void *data )
{
    if ( !H5Dexists( fid, name ) ) {
        AMP_ERROR( "Variable does not exist in file: " + std::string( name ) );
        return;
    }
    hid_t dataset   = H5Dopen2( fid, name.data(), H5P_DEFAULT );
    hid_t dataspace = H5Dget_space( dataset );
    hsize_t dims0[10];
    int ndim = H5Sget_simple_extent_dims( dataspace, dims0, NULL );
    AMP_ASSERT( ndim == 1 && dims0[0] == N_bytes );
    H5Dread( dataset, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
    H5Dclose( dataset );
    H5Sclose( dataspace );
}


#else // No HDF5
// Dummy implementations for no HDF5
hid_t openHDF5( const std::string_view &, const char *, Compression ) { return 0; }
void closeHDF5( hid_t ) {}
bool H5Gexists( hid_t, const std::string_view & ) { return false; }
bool H5Dexists( hid_t, const std::string_view & ) { return false; }
hid_t createGroup( hid_t, const std::string_view & ) { return 0; }
hid_t openGroup( hid_t, const std::string_view & ) { return 0; }
void closeGroup( hid_t ) {}
void writeHDF5( hid_t, const std::string_view &, size_t, const std::bytes * );
void readHDF5( hid_t, const std::string_view &, size_t, std::bytes * );
#endif

} // namespace AMP
