#include "AMP/IO/FileSystem.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UtilityMacros.h"

#include <algorithm>
#include <fstream>


/****************************************************************************
 *  Filesystem utilities                                                     *
 ****************************************************************************/
bool AMP::IO::fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}

std::string AMP::IO::path( const std::string &filename )
{
    size_t pos = 0;
    if ( filename.find_last_of( 47 ) != std::string::npos )
        pos = filename.find_last_of( 47 );
    if ( filename.find_last_of( 92 ) != std::string::npos )
        pos = std::max( pos, filename.find_last_of( 92 ) );
    return filename.substr( 0, pos );
}

std::string AMP::IO::filename( const std::string &filename )
{
    size_t pos = 0;
    if ( filename.find_last_of( 47 ) != std::string::npos )
        pos = filename.find_last_of( 47 );
    if ( filename.find_last_of( 92 ) != std::string::npos )
        pos = std::max( pos, filename.find_last_of( 92 ) );
    if ( pos != 0 )
        return filename.substr( pos + 1 );
    return filename;
}

void AMP::IO::renameFile( const std::string &old_filename, const std::string &new_filename )
{
    AMP_ASSERT( !old_filename.empty() );
    AMP_ASSERT( !new_filename.empty() );
    rename( old_filename.c_str(), new_filename.c_str() );
}

void AMP::IO::deleteFile( const std::string &filename )
{
    AMP_ASSERT( !filename.empty() );
    if ( fileExists( filename ) ) {
        int error = remove( filename.c_str() );
        AMP_INSIST( error == 0, "Error deleting file" );
    }
}
size_t AMP::IO::fileSize( const std::string &filename )
{
    AMP_ASSERT( !filename.empty() );
    if ( !fileExists( filename ) )
        return 0;
    auto f = fopen( filename.data(), "rb" );
    fseek( f, 0, SEEK_END );
    size_t bytes = ftell( f );
    fclose( f );
    return bytes;
}
std::string AMP::IO::getSuffix( const std::string &filename )
{
    size_t pos = filename.rfind( '.' );
    if ( pos == std::string::npos )
        return std::string();
    std::string suffix = filename.substr( pos + 1 );
    for ( auto &c : suffix )
        c = std::tolower( c );
    return suffix;
}

void AMP::IO::recursiveMkdir( const std::string &path, mode_t mode, bool only_node_zero_creates )
{
    AMP_MPI comm = AMP_MPI( AMP_COMM_WORLD );
    // Create/check a directory
    auto check = []( const std::string &path ) {
        struct stat status;
        if ( stat( path.c_str(), &status ) == 0 ) {
            if ( S_ISDIR( status.st_mode ) ) {
                return true;
            } else {
                AMP_ERROR( "Cannot create directory \"" + path +
                           "\" because it exists and is NOT a directory" );
            }
        }
        return false;
    };
    auto create = [check, mode]( const std::string &path ) {
        if ( check( path ) )
            return;
        if ( mkdir( path.c_str(), mode ) != 0 )
            AMP_ERROR( "Cannot create directory " + path );
    };
    // Create the parent directories (if they exist)
    bool createPath = comm.getRank() == 0 || !only_node_zero_creates;
    if ( !path.empty() && createPath ) {
        auto is_dir = []( char c ) { return c == '/' || c == 92; };
        auto it     = std::find_if( path.begin(), path.end(), is_dir );
        for ( ; it != path.end(); ++it ) {
            auto path2 = path.substr( 0, std::distance( path.begin(), it ) );
            create( path2 );
        }
        create( path );
    }
    // Make sure all processors wait until node zero creates the directory structure.
    if ( only_node_zero_creates )
        comm.barrier();
}
