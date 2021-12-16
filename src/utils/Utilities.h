#ifndef included_AMP_Utilities
#define included_AMP_Utilities


#include "AMP/utils/UtilityMacros.h"

#include "StackTrace/Utilities.h"

#include <chrono>
#include <cstdarg>
#include <limits>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <string>
#include <string_view>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>


namespace AMP {


class Database;

#ifdef _MSC_VER
#include <direct.h>
typedef int mode_t;
#define S_ISDIR( m ) ( ( (m) &S_IFMT ) == S_IFDIR )
#define S_IRUSR 0
#define S_IWUSR 0
#define S_IXUSR 0
#endif


// \cond HIDDEN_SYMBOLS
template<class T>
inline T type_default_tol()
{
    if constexpr ( std::is_integral<T>::value )
        return 0;
    else if constexpr ( std::is_same<T, double>::value )
        return 1e-12;
    else if constexpr ( std::is_floating_point<T>::value )
        return pow( std::numeric_limits<T>::epsilon(), (T) 0.77 );
    else
        return T();
}
// \endcond


/*!
 * Utilities is a namespace containing basic routines for error
 * reporting, file manipulations, etc.  Included are a set of \ref Macros "macros" that are commonly
 * used.
 */
namespace Utilities {


// Include functions from StackTrace
using StackTrace::Utilities::abort;
using StackTrace::Utilities::exec;
using StackTrace::Utilities::getMemoryUsage;
using StackTrace::Utilities::getSystemMemory;
using StackTrace::Utilities::tick;
using StackTrace::Utilities::time;


//! Enum for the operating system
enum class OS { macOS, Linux, Windows, Unknown };

//! Return the OS
constexpr OS getOS();


/*!
 * Set an environmental variable
 * @param name              The name of the environmental variable
 * @param value             The value to set
 */
void setenv( const char *name, const char *value );

/*!
 * Get an environmental variable
 * @param name              The name of the environmental variable
 * @return                  The value of the enviornmental variable
 */
std::string getenv( const char *name );

/*
 * Create the directory specified by the path string.  Permissions are set
 * by default to rwx by user.  The intermediate directories in the
 * path are created if they do not already exist.  When
 * only_node_zero_creates is true, only node zero creates the
 * directories.  Otherwise, all nodes create the directories.
 */
void recursiveMkdir( const std::string &path,
                     mode_t mode                 = ( S_IRUSR | S_IWUSR | S_IXUSR ),
                     bool only_node_zero_creates = true );


//! Return the path to the file
std::string path( const std::string &filename );

//! Return the filename (strip the path)
std::string filename( const std::string &filename );


//! Check if a file exists and return true if it does
bool fileExists( const std::string &filename );

//! Return the file size
size_t fileSize( const std::string &filename );


//! Rename a file from old file name to new file name.
void renameFile( const std::string &old_filename, const std::string &new_filename );


//! Delete a file.  If the file does not exist, nothing will happen.
void deleteFile( const std::string &filename );


//! Get the lower case suffix for a file
std::string getSuffix( const std::string &filename );


/*!
 * Convert an integer to a string.
 *
 * The returned string is padded with zeros as needed so that it
 * contains at least the number of characters indicated by the
 * minimum width argument.  When the number is positive, the
 * string is padded on the left. When the number is negative,
 * the '-' sign appears first, followed by the integer value
 * padded on the left with zeros.  For example, the statement
 * intToString(12, 5) returns "00012" and the statement
 * intToString(-12, 5) returns "-0012".
 */
std::string intToString( int num, int min_width = 1 );


/*!
 * Convert common integer values to strings.
 *
 * These are simply wrappers around intToString that ensure the
 * same width is uniformally used when converting to string
 * representations.
 */
std::string nodeToString( int num );
std::string processorToString( int num );
std::string patchToString( int num );
std::string levelToString( int num );
std::string blockToString( int num );


/*!
 * Soft equal checks if two numbers are within the given precision
 * True iff abs(v1-v2)/v1 < tol
 * \param v1     scalar floating point value
 * \param v2     scalar floating point value
 * \param tol    relative tolerance
 */
template<class T>
inline bool approx_equal( const T &v1, const T &v2, const T tol = type_default_tol<T>() )
{
    // Compute the absolute tolerance
    T tol2 = static_cast<T>( tol * std::max( fabs( (double) ( v1 ) ), fabs( (double) ( v2 ) ) ) );
    // Check if the two value are less than tolerance
    return fabs( (double) ( v1 - v2 ) ) <= tol2;
}

/*!
 * Soft equal checks if two numbers are equivalent within the given precision
 * True iff abs(v1-v2) < tol
 * \param v1     scalar floating point value
 * \param v2     scalar floating point value
 * \param tol    relative tolerance
 */
template<class T>
inline bool approx_equal_abs( const T &v1, const T &v2, const T tol = type_default_tol<T>() )
{
    return fabs( (double) ( v1 - v2 ) ) <= tol; // Check if the two value are less than tolerance
}


/*!
 * Quicksort a std::vector
 * \param N      Number of entries to sort
 * \param x      vector to sort
 */
template<class T>
void quicksort( size_t N, T *x );

/*!
 * Quicksort a std::vector
 * \param x      vector to sort
 */
template<class T>
inline void quicksort( std::vector<T> &x )
{
    quicksort( x.size(), x.data() );
}

/*!
 * Quicksort a std::vector
 * \param N      Number of entries to sort
 * \param x      Vector to sort
 * \param y      Extra values to be sorted with X
 */
template<class T1, class T2>
void quicksort( size_t N, T1 *x, T2 *y );

/*!
 * Quicksort a std::vector
 * \param x      Vector to sort
 * \param y      Extra values to be sorted with X
 */
template<class T1, class T2>
inline void quicksort( std::vector<T1> &x, std::vector<T2> &y )
{
    if ( x.size() != y.size() )
        AMP_ERROR( "x and y must be the same size" );
    quicksort( x.size(), x.data(), y.data() );
}


/*!
 * Get the unique set on a std::vector
 * \param x      vector to create the unique set (elements will be returned in sorted order)
 */
template<class T>
void unique( std::vector<T> &x );

/*!
 * Subroutine to perform the unique operation on the elements in X
 * This function performs the unique operation on the values in X storing them in Y.
 * It also returns the index vectors I and J such that Y[k] = X[I[k]] and X[k] = Y[J[k]].
 * @param X         Points to sort (nx)
 * @param I         The index vector I (ny)
 * @param J         The index vector J (nx)
 */
template<class T>
void unique( std::vector<T> &X, std::vector<size_t> &I, std::vector<size_t> &J );


/*!
 * Search a std::vector for the first entry >= the given value
 * This routine only works on sorted arrays and does not check if the array is sorted
 * This routine returns the size of the vector if no entries in the vector are >= the desired entry.
 * \param N      Number of entires to search
 * \param x      vector to sort
 * \param value  Value to search for
 */
template<class T>
size_t findfirst( size_t N, const T *x, const T &value );

/*!
 * Search a std::vector for the first entry >= the given value
 * This routine only works on sorted arrays and does not check if the array is sorted
 * This routine returns the size of the vector if no entries in the vector are >= the desired entry.
 * \param x      vector to sort
 * \param value  Value to search for
 */
template<class T>
inline size_t findfirst( const std::vector<T> &x, const T &value )
{
    return findfirst( x.size(), x.data(), value );
}


/*!
 * Function to perform linear interpolation
 * \param x     x-coordinates
 * \param f     function values at the coordinates ( Nx )
 * \param xi    x-coordinate of desired point
 */
double linear( const std::vector<double> &x, const std::vector<double> &f, double xi );

/*!
 * Function to perform tri-linear interpolation
 * \param x     x-coordinates
 * \param y     y-coordinates
 * \param f     function values at the coordinates ( Nx x Ny )
 * \param xi    x-coordinate of desired point
 * \param yi    y-coordinate of desired point
 */
double bilinear( const std::vector<double> &x,
                 const std::vector<double> &y,
                 const std::vector<double> &f,
                 double xi,
                 double yi );

/*!
 * Function to perform tri-linear interpolation
 * \param x     x-coordinates
 * \param y     y-coordinates
 * \param z     z-coordinates
 * \param f     function values at the coordinates ( Nx x Ny x Nz )
 * \param xi    x-coordinate of desired point
 * \param yi    y-coordinate of desired point
 * \param zi    z-coordinate of desired point
 */
double trilinear( const std::vector<double> &x,
                  const std::vector<double> &y,
                  const std::vector<double> &z,
                  const std::vector<double> &f,
                  double xi,
                  double yi,
                  double zi );

//! Create a hash key from a char array
constexpr unsigned int hash_char( const char * );


// Function to demangle a string (e.g. from typeid)
std::string demangle( const std::string &name );


//! Get the type name (does not match typeid, does not work for all compilers)
template<typename T>
constexpr std::string_view type_name();

//! Get the type hash (does not match typeid, does not work for all compilers)
template<typename T>
constexpr unsigned int type_hash();


//! Get the prime factors for a number
std::vector<int> factor( uint64_t );


//! Check if a number is prime
bool isPrime( uint64_t );


//! Return all prime numbers <= x
std::vector<uint64_t> primes( uint64_t );


/*!
 * Sleep for X ms
 * @param N         Time to sleep (ms)
 */
inline void sleep_ms( int N ) { std::this_thread::sleep_for( std::chrono::milliseconds( N ) ); }

/*!
 * Sleep for X s
 * @param N         Time to sleep (s)
 */
inline void sleep_s( int N ) { std::this_thread::sleep_for( std::chrono::seconds( N ) ); }

//! Print AMP Banner
void printBanner();

//! Null use function
void nullUse( void * );

//! std::string version of sprintf
inline std::string stringf( const char *format, ... );

//! Print a vector
template<class TYPE>
std::string to_string( const std::vector<TYPE> &x );


//! Print a database
[[deprecated( "This function will be removed soon, use Database::print" )]] void
printDatabase( const Database &, std::ostream &, const std::string &indent = "" );


//! Stack based vector
template<class TYPE, std::size_t CAPACITY>
class stackVector final
{
public:
    stackVector() : d_size( 0 ) {}
    size_t size() const { return d_size; }
    bool empty() const { return d_size == 0; }
    void push_back( const TYPE &v )
    {
        if ( d_size == CAPACITY )
            throw std::bad_alloc();
        d_data[d_size++] = v;
    }
    void clear()
    {
        for ( size_t i = 0; i < d_size; i++ )
            d_data[i] = TYPE();
        d_size = 0;
    }
    TYPE &operator[]( size_t i ) { return d_data[i]; }
    TYPE *begin() { return d_data; }
    TYPE *end() { return d_data + d_size; }
    TYPE &back() { return d_data[d_size - 1]; }
    TYPE *data() { return d_size == 0 ? nullptr : d_data; }
    void pop_back() { d_size = std::max<size_t>( d_size, 1 ) - 1; }
    const TYPE &operator[]( size_t i ) const { return d_data[i]; }
    const TYPE *begin() const { return d_data; }
    const TYPE *end() const { return d_data + d_size; }
    const TYPE &back() const { return d_data[d_size - 1]; }
    template<class... Args>
    void emplace_back( Args &&...args )
    {
        if ( d_size == CAPACITY )
            throw std::bad_alloc();
        d_data[d_size++] = TYPE( args... );
    }

private:
    uint32_t d_size;
    TYPE d_data[CAPACITY];
};


} // namespace Utilities
} // namespace AMP


#include "AMP/utils/Utilities.hpp"


#endif
