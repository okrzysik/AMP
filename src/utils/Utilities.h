#ifndef included_AMP_Utilities
#define included_AMP_Utilities


#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS		// Supress depreciated warnings for visual studio
#endif


// Include the utility macros
#include "UtilityMacros.h"


#include <string>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "Logger.h"
#include <typeinfo>
#include <limits>

namespace AMP {
  

#ifdef _MSC_VER
    #define _CRT_SECURE_NO_WARNINGS
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <direct.h>
    typedef int mode_t;
    #define  S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
    #define S_IRUSR 0
    #define S_IWUSR 0
    #define S_IXUSR 0
#endif


// \cond HIDDEN_SYMBOLS
template<class T>  inline T type_default_tol();
template<>  inline int type_default_tol<int>() { return 0; }
template<>  inline unsigned int type_default_tol<unsigned int>() { return 0; }
template<>  inline size_t type_default_tol<size_t>() { return 0; }
template<>  inline double type_default_tol<double>() { return 1e-12; }
template<class T>  inline T type_default_tol() { return pow(std::numeric_limits<T>::epsilon(),(T)0.77); }
// \endcond





/*!
 * Utilities is a Singleton class containing basic routines for error 
 * reporting, file manipulations, etc.  Included are a set of \ref Macros "macros" that are commonly used.
 */
namespace Utilities
{
    /*!
     * Create the directory specified by the path string.  Permissions are set 
     * by default to rwx by user.  The intermediate directories in the 
     * path are created if they do not already exist.  When 
     * only_node_zero_creates is true, only node zero creates the 
     * directories.  Otherwise, all nodes create the directories.
     */
    void recursiveMkdir(const std::string& path, 
			      mode_t mode = (S_IRUSR|S_IWUSR|S_IXUSR),
			      bool only_node_zero_creates = true);

    /*!
     * Check if a file exists and return true if it does
     */
    bool fileExists( const std::string& filename );


    /*!
     * Rename a file from old file name to new file name.
     */
    void renameFile(const std::string& old_filename, 
                          const std::string& new_filename);


    /*!
     * Delete a file.  If the file does not exist, nothing will happen.
     */
    void deleteFile( const std::string& filename );


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
    std::string intToString(int num, int min_width = 1);

    /*!
     * Convert common integer values to strings.
     *
     * These are simply wrappers around intToString that ensure the 
     * same width is uniformally used when converting to string
     * representations.
     */
    std::string nodeToString(int num);
    std::string processorToString(int num);
    std::string patchToString(int num);
    std::string levelToString(int num);
    std::string blockToString(int num);

    /*!
     * Aborts the run after printing an error message with file and
     * linenumber information.
     */
    void abort(const std::string &message, 
		     const std::string &filename,
		     const int line);

    /*!
     * Soft equal checks if two numbers are within the given precision
     * True iff abs(v1-v2)/v1 < tol
     * \param v1     scalar floating point value
     * \param v2     scalar floating point value
     * \param tol    relative tolerance
     */
    template<class T>
    inline bool approx_equal(const T &v1, const T &v2, const T tol = type_default_tol<T>() ) {
        T tol2 = tol*std::max( fabs(v1), fabs(v2) );    // Compute the absolute tolerance
        return fabs(v1-v2)<=tol2;                       // Check if the two value are less than tolerance
    }

    /*!
     * Soft equal checks if two numbers are equivalent within the given precision
     * True iff abs(v1-v2) < tol
     * \param v1     scalar floating point value
     * \param v2     scalar floating point value
     * \param tol    relative tolerance
     */
    template<class T>
    inline bool approx_equal_abs(const T &v1, const T &v2, const T tol = type_default_tol<T>() ) {
        return fabs(v1-v2)<=tol;                    // Check if the two value are less than tolerance
    }

    /*!
     * Quicksort a std::vector
     * \param x      vector to sort
     */
    template<class T>
    void quicksort(std::vector<T> &x);

    /*!
     * Quicksort a std::vector
     * \param x      Vector to sort
     * \param y      Extra values to be sorted with X
     */
    template<class T1, class T2>
    void quicksort(std::vector<T1> &x, std::vector<T2> &y);

    /*!
     * Get the unique set on a std::vector
     * \param x      vector to create the unique set (elements will be returned in sorted order)
     */
    template<class T>
    void unique(std::vector<T> &x);

    /*!
     * Search a std::vector for the first entry >= the given value
     * This routine only works on sorted arrays and does not check if the array is sorted
     * This routine returns the size of the vector if no entries in the vector are >= the desired entry.
     * \param x      vector to sort
     * \param value  Value to search for
     */
    template<class T>
    size_t findfirst(const std::vector<T> &x, const T &value);

    /*!
     * Function to perform linear interpolation
     * \param x     x-coordinates
     * \param f     function values at the coordinates ( Nx )
     * \param xi    x-coordinate of desired point
     */
    double linear( const std::vector<double>& x, const std::vector<double>& f, double xi );

    /*!
     * Function to perform tri-linear interpolation
     * \param x     x-coordinates
     * \param y     y-coordinates
     * \param f     function values at the coordinates ( Nx x Ny )
     * \param xi    x-coordinate of desired point
     * \param yi    y-coordinate of desired point
     */
    double bilinear( const std::vector<double>& x, const std::vector<double>& y, 
        const std::vector<double>& f, double xi, double yi );

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
    double trilinear( const std::vector<double>& x, const std::vector<double>& y, 
        const std::vector<double>& z, const std::vector<double>& f, double xi, double yi, double zi );

    //! Create a hash key from a char array
    unsigned int hash_char(const char*);

    //! Get the prime factors for a number
    std::vector<int> factor(size_t);

    /*!
     * Function to get the memory usage.
     * This function will return the total memory used by the application.
     * Note: depending on the implimentation, this number may be rounded to
     * to a multiple of 1024 (2^3n).
     * If this function fails, it will return 0.
     */
    size_t getMemoryUsage();

    //! Function to get the current call stack
    std::vector<std::string> getCallStack();

    //! Print AMP Banner
    void printBanner();


    //! Triplet version of std::pair
    template<class A, class B, class C>
    struct triplet{ 
        A first; 
        B second; 
        C third; 
        inline triplet() {
            first = A();
            second = B();
            third = C();
        }
        inline triplet(A a, B b, C c) {
            first = a;
            second = b;
            third = c;
        }
        inline triplet(const triplet& rhs ) {
            first = rhs.first;
            second = rhs.second;
            third = rhs.third;
        }
        inline bool operator== (const triplet& rhs ) const {
            return first==rhs.first && second==rhs.second && third==rhs.third;
        }
        inline bool operator!= (const triplet& rhs ) const {
            return first!=rhs.first || second!=rhs.second || third!=rhs.third;
        }
        inline bool operator>= (const triplet& rhs ) const {
            if ( first < rhs.first )        { return false; }
            else if ( first > rhs.first )   { return true;  }
            if ( second < rhs.second )      { return false; }
            else if ( second > rhs.second ) { return true;  }
            return third>=rhs.third;
        }
        inline bool operator> (const triplet& rhs ) const {
            if ( first < rhs.first )        { return false; }
            else if ( first > rhs.first )   { return true;  }
            if ( second < rhs.second )      { return false; }
            else if ( second > rhs.second ) { return true;  }
            return third>rhs.third;
        }
        inline bool operator< (const triplet& rhs ) const {
            if ( first > rhs.first )        { return false; }
            else if ( first < rhs.first )   { return true;  }
            if ( second > rhs.second )      { return false; }
            else if ( second < rhs.second ) { return true;  }
            return third<rhs.third;
        }
        inline bool operator<= (const triplet& rhs ) const {
            if ( first > rhs.first )        { return false; }
            else if ( first < rhs.first )   { return true;  }
            if ( second > rhs.second )      { return false; }
            else if ( second < rhs.second ) { return true;  }
            return third<=rhs.third;
        }
    };

}


// templated quicksort routine
template <class T>
void Utilities::quicksort(std::vector<T> &x)
{
    int n = (int) x.size();
    if ( n <= 1 )
        return;
    T *arr = &x[0];
    bool test;
    int i, ir, j, jstack, k, l, istack[100];
    T a, tmp_a;
    jstack = 0;
    l = 0;
    ir = n-1;
    while (1) {
        if ( ir-l < 7 ) {             // Insertion sort when subarray small enough.
            for ( j=l+1; j<=ir; j++ ) {
                a = arr[j];
                test = true;
                for (i=j-1; i>=0; i--) {
                    if ( arr[i] < a ) {
                        arr[i+1] = a;
                        test = false;
                        break;
                    }
                    arr[i+1] = arr[i];
                }
                if ( test ) {
                    i = l-1;
                    arr[i+1] = a;
                }
            }
            if ( jstack==0 )
                return;
            ir = istack[jstack];    // Pop stack and begin a new round of partitioning.
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir)/2;           // Choose median of left, center and right elements as partitioning
                                    // element a. Also rearrange so that a(l) < a(l+1) < a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l+1];
            arr[l+1] = tmp_a;
            if ( arr[l]>arr[ir] ) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l+1] > arr[ir] ) {
                tmp_a = arr[l+1];
                arr[l+1] = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l] > arr[l+1] ) {
                tmp_a = arr[l];
                arr[l] = arr[l+1];
                arr[l+1] = tmp_a;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l+1];           // Partitioning element.
            for (i=l+2; i<=ir; i++) { 
                if ( arr[i]<a ) 
                    continue;
                while ( arr[j]>a )  // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;          // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i];     // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
            }
            arr[l+1] = arr[j];      // Insert partitioning element in both arrays.
            arr[j] = a;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir-i+1 >= j-l ) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}


// templated quicksort routine
template <class T1, class T2>
void Utilities::quicksort(std::vector<T1> &x, std::vector<T2> &y)
{
    if ( x.size() != y.size() )
        AMP_ERROR("x and y must be the same size");
    int n = (int) x.size();
    if ( n <= 1 )
        return;
    T1 *arr = &x[0];
    T2 *brr = &y[0];
    bool test;
    int i, ir, j, jstack, k, l, istack[100];
    T1 a, tmp_a;
    T2 b, tmp_b;
    jstack = 0;
    l = 0;
    ir = n-1;
    while (1) {
        if ( ir-l < 7 ) {             // Insertion sort when subarray small enough.
            for ( j=l+1; j<=ir; j++ ) {
                a = arr[j];
                b = brr[j];
                test = true;
                for (i=j-1; i>=0; i--) {
                    if ( arr[i] < a ) {
                        arr[i+1] = a;
                        brr[i+1] = b;
                        test = false;
                        break;
                    }
                    arr[i+1] = arr[i];
                    brr[i+1] = brr[i];
                }
                if ( test ) {
                    i = l-1;
                    arr[i+1] = a;
                    brr[i+1] = b;
                }
            }
            if ( jstack==0 )
                return;
            ir = istack[jstack];    // Pop stack and begin a new round of partitioning.
            l = istack[jstack-1];
            jstack -= 2;
        } else {
            k = (l+ir)/2;           // Choose median of left, center and right elements as partitioning
                                    // element a. Also rearrange so that a(l) ? a(l+1) ? a(ir).
            tmp_a = arr[k];
            arr[k] = arr[l+1];
            arr[l+1] = tmp_a;
            tmp_b = brr[k];
            brr[k] = brr[l+1];
            brr[l+1] = tmp_b;
            if ( arr[l]>arr[ir] ) {
                tmp_a = arr[l];
                arr[l] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[ir];
                brr[ir] = tmp_b;
            }
            if ( arr[l+1] > arr[ir] ) {
                tmp_a = arr[l+1];
                arr[l+1] = arr[ir];
                arr[ir] = tmp_a;
                tmp_b = brr[l+1];
                brr[l+1] = brr[ir];
                brr[ir] = tmp_b;
            }
            if ( arr[l] > arr[l+1] ) {
                tmp_a = arr[l];
                arr[l] = arr[l+1];
                arr[l+1] = tmp_a;
                tmp_b = brr[l];
                brr[l] = brr[l+1];
                brr[l+1] = tmp_b;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l+1];           // Partitioning element.
            b = brr[l+1];
            for (i=l+2; i<=ir; i++) { 
                if ( arr[i]<a ) 
                    continue;
                while ( arr[j]>a )  // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;          // Pointers crossed. Exit with partitioning complete.
                tmp_a = arr[i];     // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
                tmp_b = brr[i];
                brr[i] = brr[j];
                brr[j] = tmp_b;
            }
            arr[l+1] = arr[j];      // Insert partitioning element in both arrays.
            arr[j] = a;
            brr[l+1] = brr[j];
            brr[j] = b;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir-i+1 >= j-l ) {
                istack[jstack] = ir;
                istack[jstack-1] = i;
                ir = j-1;
            } else {
                istack[jstack] = j-1;
                istack[jstack-1] = l;
                l = i;
            }
        }
    }
}


/************************************************************************
* Subroutine to find the unique elements in a list                      *
************************************************************************/
template <class T>
void Utilities::unique(std::vector<T> &x)
{
    if ( x.size()==0 )
        return;
    // First perform a quicksort
    Utilities::quicksort(x);
    // Next remove duplicate entries
    size_t pos = 1;
    for (size_t i=1; i<x.size(); i++) {
        if ( x[i] != x[pos-1] ) {
            x[pos] = x[i];
            pos++;
        }
    }
    if ( pos < x.size() )
        x.resize(pos);
}


/************************************************************************
* Subroutine to find the first element in X which is greater than Y     *
* using a simple hashing technique.  This is the a faster method, but   *
* requires the vector X to be in ascending order.                       *
* Returns -1 if no value is larger.                                     *
************************************************************************/
template<class T>
size_t Utilities::findfirst(const std::vector<T> &x_in, const T &value) 
{
    size_t n = x_in.size();
    AMP_INSIST(n>0,"x must not be empty");
    const T *x = &x_in[0];   // Use the pointer for speed
    // Check if value is within the range of x
    if ( value <= x[0] )
        return 0;
    else if ( value > x[n-1] )
        return n;
    // Perform the search
    size_t lower = 0;
    size_t upper = n-1;
    size_t index;
    while ( (upper-lower) != 1 ) {
        index = (upper+lower)/2;
        if ( x[index] >= value )
            upper = index;
        else
            lower = index;
    }
    index = upper;
    return index;
}


}

#endif


