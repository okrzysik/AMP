// This file contains helper functions and interfaces for filesystem access
// These should be replaced with std::filesystem in the future
#ifndef included_AMP_FileSystem
#define included_AMP_FileSystem


#include <stdio.h>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <sys/types.h>


#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 ) || \
    defined( _MSC_VER )
    #include <direct.h>
typedef int mode_t;
    #define S_ISDIR( m ) ( ( (m) &S_IFMT ) == S_IFDIR )
    #define S_IRUSR 0
    #define S_IWUSR 0
    #define S_IXUSR 0
    #define S_IRWXU 0
#endif


namespace AMP::IO {


/*
 * Create the directory specified by the path string.  Permissions are set
 * by default to rwx by user.  The intermediate directories in the
 * path are created if they do not already exist.
 */
void recursiveMkdir( const std::string &path, mode_t mode = ( S_IRUSR | S_IWUSR | S_IXUSR ) );
[[deprecated( "No need to specify only_node_zero_creates, use recursiveMkdir(path,mode)" )]] void
recursiveMkdir( const std::string &, mode_t, bool );


//! Return the path to the file
std::string path( const std::string &filename );

//! Return the filename (strip the path)
std::string filename( const std::string &filename );


//! Check if a file exists and return true if it does
bool fileExists( const std::string &filename );

//! Return the file size
size_t fileSize( const std::string &filename );


//! Rename a file from old file name to new file name
void renameFile( const std::string &old_filename, const std::string &new_filename );


//! Delete a file.  If the file does not exist, nothing will happen
void deleteFile( const std::string &filename );


//! Get the lower case suffix for a file
std::string getSuffix( const std::string &filename );


} // namespace AMP::IO

#endif
