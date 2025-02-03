#include "AMP/IO/FileSystem.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <thread>


#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || \
    defined( _MSC_VER )
int mkdir( const char *path, mode_t ) { return _mkdir( path ); }
#endif


/****************************************************************************
 *  Filesystem utilities                                                     *
 ****************************************************************************/
constexpr static char slash[] = { 0x2F, 0x5C, 0x0 };
bool AMP::IO::fileExists( const std::string &filename )
{
    std::ifstream ifile( filename.c_str() );
    return ifile.good();
}

std::string AMP::IO::path( const std::string &filename )
{
    size_t pos = filename.find_last_of( slash );
    if ( pos == std::string::npos )
        return {};
    return filename.substr( 0, pos );
}

std::string AMP::IO::filename( const std::string &filename )
{
    size_t pos = filename.find_last_of( slash );
    if ( pos == std::string::npos )
        return filename;
    return filename.substr( pos + 1 );
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
void AMP::IO::recursiveMkdir( const std::string &path, mode_t mode, bool )
{
    recursiveMkdir( path, mode );
}
void AMP::IO::recursiveMkdir( const std::string &path, mode_t mode )
{
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
        if ( mkdir( path.c_str(), mode ) != 0 ) {
            if ( errno == EEXIST && check( path ) ) {
                // Another rank probably created the directory first
                return;
            }
            AMP_ERROR( "Cannot create directory " + path + " (" +
                       AMP::Utilities::getLastErrnoString() + ")" );
        }
    };
    // Create the parent directories (if they exist)
    if ( path.empty() )
        return;
    if ( AMP_MPI( AMP_COMM_WORLD ).getRank() != 0 )
        std::this_thread::sleep_for( std::chrono::milliseconds( 20 ) );
    std::string split = "\\/";
    size_t pos        = path.find_first_of( split, 0 );
    while ( pos != std::string::npos ) {
        create( path.substr( 0, pos ) );
        pos = path.find_first_of( split, pos + 1 );
    }
    create( path );
}
