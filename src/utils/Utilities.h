//
// File:        $URL: file:///usr/casc/samrai/repository/AMP/tags/v-2-4-4/source/toolbox/base/Utilities.h $
// Package:     AMP toolbox
// Copyright:   (c) 1997-2008 Lawrence Livermore National Security, LLC
// Revision:    $LastChangedRevision: 2249 $
// Modified:    $LastChangedDate: 2008-07-03 08:17:20 -0700 (Thu, 03 Jul 2008) $
// Description: Utility functions for error reporting, file manipulation, etc.
//

#ifndef included_AMP_Utilities
#define included_AMP_Utilities


#include <string>
#include <math.h>
#include <vector>
#include "IOStream.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "Logger.h"
#include <typeinfo>

namespace AMP {
  
#ifdef _MSC_VER
    #include <sys/types.h>
    #include <sys/stat.h>
    #include <direct.h>
    typedef int mode_t;
    #define  S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
    #define S_IRUSR 0
    #define S_IWUSR 0
    #define S_IXUSR 0
#endif

/*!
 * Utilities is a Singleton class containing basic routines for error 
 * reporting, file manipulations, etc.
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
     * Rename a file from old file name to new file name.
     */
    void renameFile(const std::string& old_filename, 
                          const std::string& new_filename);

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
     * \param v1     scalar floating point value
     * \param v2     scalar floating point value
     * \param tol    relative tolerance
     */
    template<class T>
    inline bool approx_equal(const T &v1, const T &v2, const T tol=1.0e-12) {
        T tol2 = tol*fabs(v2);                      // Compute the absolute tolerance
        if ( tol2 < 1.0e-14 ) { tol2 = 1.0e-14; }   // Fix the minimum tolerance
        return fabs(v1-v2)<tol2;                    // Check if the two value are less than tolerance
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
     * Search a std::vector for the first entry >= the given value
     * This routine only works on sorted arrays and does not check if the array is sorted
     * This routine returns the size of the vector if no entries in the vector are >= the desired entry.
     * \param x      vector to sort
     * \param value  Value to search for
     */
    template<class T>
    size_t findfirst(const std::vector<T> &x, const T &value);

    //! Create a hash key from a char array
    unsigned int hash_char(const char*);

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

}


/*!
 * A statement that does nothing, for insure++ make it something 
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if(0) int nullstatement=0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a variable, use to avoid GNU compiler 
 * warnings about unused variables.
 */
#define NULL_USE(variable) do { \
       if(0) {char *temp = (char *)&variable; temp++;} \
    } while (0)

/*!
 * Throw an error exception from within any C++ source code.  The 
 * macro argument may be any standard ostream expression.  The file and
 * line number of the abort are also printed.
 */
#ifndef LACKS_SSTREAM
#define AMP_ERROR(X) do {					\
      std::ostringstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define AMP_ERROR(X) do {					\
      std::ostrstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /*!
    * Print a warning without exit.  Print file and line number of the warning.
    */
#ifndef LACKS_SSTREAM
#define AMP_WARNING(X) do {					\
      std::ostringstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Logger::getInstance() -> logWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define AMP_WARNING(X) do {					\
      std::ostrstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Logger::getInstance() -> logWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif


/*!
 * Print a debug without exit.  Print file and line number of the debug.
 * \todo 
 * Fix AMP_DEBUG so it uses a consistent stream io with AMP
 */
#ifndef LACKS_SSTREAM
#define AMP_DEBUG(X) do {					\
      std::ostringstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Logger::getInstance() -> logDebug(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define AMP_DEBUG(X) do {					\
      std::ostrstream tboxos;					\
      tboxos << X << std::ends;					\
      AMP::Logger::getInstance() -> logDebug(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif


/*! \def AMP_ASSERT
 * \brief Assert error
 * Throw an error exception from within any C++ source code if the
 * given expression is not true.  This is a parallel-friendly version
 * of assert.
 * The file and line number of the abort are also printed.
 * \todo 
 * It would be helpful to add a full stack trace
 */
/*! \def AMP_INSIST
 * \brief Insist error
 * Throw an error exception from within any C++ source code if the
 * given expression is not true.  This will also print the given message.
 * This is a parallel-friendly version of assert.
 * The file and line number of the abort are also printed.
 * \todo 
 * It would be helpful to add a full stack trace
 */
#ifdef HAVE_STRINGIZE
    #ifndef LACKS_SSTREAM

        #define AMP_ASSERT(EXP) do {                                       \
            if ( !(EXP) ) {                                                 \
                std::ostringstream tboxos;                                  \
                tboxos << "Failed assertion: " << #EXP << std::ends;        \
                AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
            }                                                               \
        } while (0)
        #define AMP_INSIST(EXP,MSG) do {                                   \
            if ( !(EXP) ) {                                                 \
                std::ostringstream tboxos;                                  \
                tboxos << "Failed insist: " << #EXP << std::endl;           \
                tboxos << "Message: " << MSG << std::ends;                  \
                AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
            }                                                               \
        } while (0)
    #else
        #define AMP_ASSERT(EXP) do {                                       \
            if ( !(EXP) ) {                                                 \
                std::ostrstream tboxos;                                     \
                tboxos << "Failed assertion: " << #EXP << std::ends;        \
                AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
            }                                                               \
        } while (0)
        #define AMP_INSIST(EXP,MSG) do {                                   \
            if ( !(EXP) ) {                                                 \
                std::ostrstream tboxos;                                     \
                tboxos << "Failed insist: " << #EXP << std::endl;           \
                tboxos << "Message: " << MSG << std::ends;                  \
                AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);    \
            }                                                               \
        } while (0)
    #endif
#else
    #ifndef LACKS_SSTREAM
        #define AMP_ASSERT(EXP) do {                                       \
              if ( !(EXP) ) {                                               \
                 std::ostringstream tboxos;                                 \
                 tboxos << "Failed assertion: " << std::ends;               \
                 AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);   \
              }                                                             \
        } while (0)
        #define AMP_INSIST(EXP,MSG) do {                                   \
              if ( !(EXP) ) {                                               \
                 std::ostringstream tboxos;                                 \
                 tboxos << "Failed insist: " << #EXP << std::endl;          \
                 tboxos << "Message: " << MSG << std::ends;                 \
                 AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);   \
              }                                                             \
        } while (0)
    #else
        #define AMP_ASSERT(EXP) do {                                       \
              if ( !(EXP) ) {                                               \
                 std::ostrstream tboxos;                                    \
                 tboxos << "Failed assertion: " << std::ends;               \
                 AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);   \
              }                                                             \
        } while (0)
        #define AMP_INSIST(EXP,MSG) do {                                   \
              if ( !(EXP) ) {                                               \
                 std::ostrstream tboxos;                                    \
                tboxos << "Failed insist: " << #EXP << std::endl;           \
                tboxos << "Message: " << MSG << std::ends;                  \
                 AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);   \
              }                                                             \
        } while (0)
    #endif
#endif


/**
 * Macro for use when assertions are to be included
 * only when debugging.
 */
#ifdef DEBUG_CHECK_ASSERTIONS
#define AMP_CHECK_ASSERT(EXP) AMP_ASSERT(EXP)
#else
#define AMP_CHECK_ASSERT(EXP) 
#endif


/**
 * Throw an error exception from within any C++ source code.  This is
 * is similar to AMP_ERROR(), but is designed to be invoked after a
 * call to a PETSc library function.  In other words, it acts similarly
 * to the PETSc CHKERRQ(ierr) macro.
 */
#ifdef HAVE_PETSC

/*
 * In the following, "CHKERRCONTINUE(ierr);" will cause PETSc to print out
 * a stack trace that led to the error; this may be useful for debugging.
 */
 
#ifndef LACKS_SSTREAM
#define PETSC_AMP_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         std::ostringstream tboxos;							\
         AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	\
      } 									\
} while (0)
#else
#define PETSC_AMP_ERROR(ierr) do {						\
      if (ierr) {                                   				\
         std::ostrstream tboxos;							\
         CHKERRCONTINUE(ierr); 							\
         AMP::Utilities::abort(tboxos.str(), __FILE__, __LINE__);	        \
      } 									\
} while (0)
#endif
#endif


//! Get a hash key from the class type (requires the RTTI (Run-time type information) to be available.
#define TYPE_HASH(X)  AMP::Utilities::hash_char(typeid(X).name())



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
