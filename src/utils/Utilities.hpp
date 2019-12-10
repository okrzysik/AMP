#ifndef included_AMP_Utilities_hpp
#define included_AMP_Utilities_hpp


#include "AMP/utils/PIO.h"


/************************************************************************
 * Functions to hash                                                     *
 ************************************************************************/
constexpr unsigned int AMP::Utilities::hash_char( const char *name )
{
    uint32_t hash   = 5381;
    unsigned char c = 0;
    while ( ( c = *name++ ) ) {
        // hash = hash * 33 ^ c
        hash = ( ( hash << 5 ) + hash ) ^ c;
    }
    return hash;
}


/************************************************************************
 * Function to wrap printf to std::string                                *
 ************************************************************************/
inline std::string AMP::Utilities::stringf( const char *format, ... )
{
    va_list ap;
    va_start( ap, format );
    char tmp[4096];
    vsprintf( tmp, format, ap );
    va_end( ap );
    return std::string( tmp );
}


/************************************************************************
 * templated quicksort routines                                          *
 ************************************************************************/
template<class T>
void AMP::Utilities::quicksort( size_t n, T *x )
{
    if ( n <= 1 )
        return;
    T *arr = &x[0];
    bool test;
    long int i, ir, j, jstack, k, l, istack[100];
    T a, tmp_a;
    jstack = 0;
    l      = 0;
    ir     = n - 1;
    while ( 1 ) {
        if ( ir - l < 7 ) { // Insertion sort when subarray small enough.
            for ( j = l + 1; j <= ir; j++ ) {
                a    = arr[j];
                test = true;
                for ( i = j - 1; i >= 0; i-- ) {
                    if ( arr[i] < a ) {
                        arr[i + 1] = a;
                        test       = false;
                        break;
                    }
                    arr[i + 1] = arr[i];
                }
                if ( test ) {
                    i          = l - 1;
                    arr[i + 1] = a;
                }
            }
            if ( jstack == 0 )
                return;
            ir = istack[jstack]; // Pop stack and begin a new round of partitioning.
            l  = istack[jstack - 1];
            jstack -= 2;
        } else {
            k = ( l + ir ) / 2; // Choose median of left, center and right elements as partitioning
                                // element a. Also rearrange so that a(l) < a(l+1) < a(ir).
            tmp_a      = arr[k];
            arr[k]     = arr[l + 1];
            arr[l + 1] = tmp_a;
            if ( arr[l] > arr[ir] ) {
                tmp_a   = arr[l];
                arr[l]  = arr[ir];
                arr[ir] = tmp_a;
            }
            if ( arr[l + 1] > arr[ir] ) {
                tmp_a      = arr[l + 1];
                arr[l + 1] = arr[ir];
                arr[ir]    = tmp_a;
            }
            if ( arr[l] > arr[l + 1] ) {
                tmp_a      = arr[l];
                arr[l]     = arr[l + 1];
                arr[l + 1] = tmp_a;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l + 1]; // Partitioning element.
            for ( i = l + 2; i <= ir; i++ ) {
                if ( arr[i] < a )
                    continue;
                while ( arr[j] > a ) // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;       // Pointers crossed. Exit with partitioning complete.
                tmp_a  = arr[i]; // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
            }
            arr[l + 1] = arr[j]; // Insert partitioning element in both arrays.
            arr[j]     = a;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir - i + 1 >= j - l ) {
                istack[jstack]     = ir;
                istack[jstack - 1] = i;
                ir                 = j - 1;
            } else {
                istack[jstack]     = j - 1;
                istack[jstack - 1] = l;
                l                  = i;
            }
        }
    }
}
template<class T1, class T2>
void AMP::Utilities::quicksort( size_t n, T1 *x, T2 *y )
{
    if ( n <= 1 )
        return;
    T1 *arr = &x[0];
    T2 *brr = &y[0];
    bool test;
    long int i, ir, j, jstack, k, l, istack[100];
    T1 a, tmp_a;
    T2 b, tmp_b;
    jstack = 0;
    l      = 0;
    ir     = n - 1;
    while ( 1 ) {
        if ( ir - l < 7 ) { // Insertion sort when subarray small enough.
            for ( j = l + 1; j <= ir; j++ ) {
                a    = arr[j];
                b    = brr[j];
                test = true;
                for ( i = j - 1; i >= 0; i-- ) {
                    if ( arr[i] < a ) {
                        arr[i + 1] = a;
                        brr[i + 1] = b;
                        test       = false;
                        break;
                    }
                    arr[i + 1] = arr[i];
                    brr[i + 1] = brr[i];
                }
                if ( test ) {
                    i          = l - 1;
                    arr[i + 1] = a;
                    brr[i + 1] = b;
                }
            }
            if ( jstack == 0 )
                return;
            ir = istack[jstack]; // Pop stack and begin a new round of partitioning.
            l  = istack[jstack - 1];
            jstack -= 2;
        } else {
            k = ( l + ir ) / 2; // Choose median of left, center and right elements as partitioning
                                // element a. Also rearrange so that a(l) ? a(l+1) ? a(ir).
            tmp_a      = arr[k];
            arr[k]     = arr[l + 1];
            arr[l + 1] = tmp_a;
            tmp_b      = brr[k];
            brr[k]     = brr[l + 1];
            brr[l + 1] = tmp_b;
            if ( arr[l] > arr[ir] ) {
                tmp_a   = arr[l];
                arr[l]  = arr[ir];
                arr[ir] = tmp_a;
                tmp_b   = brr[l];
                brr[l]  = brr[ir];
                brr[ir] = tmp_b;
            }
            if ( arr[l + 1] > arr[ir] ) {
                tmp_a      = arr[l + 1];
                arr[l + 1] = arr[ir];
                arr[ir]    = tmp_a;
                tmp_b      = brr[l + 1];
                brr[l + 1] = brr[ir];
                brr[ir]    = tmp_b;
            }
            if ( arr[l] > arr[l + 1] ) {
                tmp_a      = arr[l];
                arr[l]     = arr[l + 1];
                arr[l + 1] = tmp_a;
                tmp_b      = brr[l];
                brr[l]     = brr[l + 1];
                brr[l + 1] = tmp_b;
            }
            // Scan up to find element > a
            j = ir;
            a = arr[l + 1]; // Partitioning element.
            b = brr[l + 1];
            for ( i = l + 2; i <= ir; i++ ) {
                if ( arr[i] < a )
                    continue;
                while ( arr[j] > a ) // Scan down to find element < a.
                    j--;
                if ( j < i )
                    break;       // Pointers crossed. Exit with partitioning complete.
                tmp_a  = arr[i]; // Exchange elements of both arrays.
                arr[i] = arr[j];
                arr[j] = tmp_a;
                tmp_b  = brr[i];
                brr[i] = brr[j];
                brr[j] = tmp_b;
            }
            arr[l + 1] = arr[j]; // Insert partitioning element in both arrays.
            arr[j]     = a;
            brr[l + 1] = brr[j];
            brr[j]     = b;
            jstack += 2;
            // Push pointers to larger subarray on stack, process smaller subarray immediately.
            if ( ir - i + 1 >= j - l ) {
                istack[jstack]     = ir;
                istack[jstack - 1] = i;
                ir                 = j - 1;
            } else {
                istack[jstack]     = j - 1;
                istack[jstack - 1] = l;
                l                  = i;
            }
        }
    }
}


/************************************************************************
 * Subroutine to find the unique elements in a list                      *
 ************************************************************************/
template<class T>
void AMP::Utilities::unique( std::vector<T> &x )
{
    if ( x.size() <= 1 )
        return;
    // First perform a quicksort
    AMP::Utilities::quicksort( x );
    // Next remove duplicate entries
    size_t pos = 1;
    for ( size_t i = 1; i < x.size(); i++ ) {
        if ( x[i] != x[pos - 1] ) {
            x[pos] = x[i];
            pos++;
        }
    }
    if ( pos < x.size() )
        x.resize( pos );
}
template<class T>
void AMP::Utilities::unique( std::vector<T> &X, std::vector<size_t> &I, std::vector<size_t> &J )
{
    const size_t neg_one = static_cast<size_t>( -1 );
    J.resize( 0 );
    J.resize( X.size(), neg_one );
    // Initialize the index vector I
    I.resize( X.size() );
    for ( size_t i = 0; i < X.size(); i++ )
        I[i] = i;
    // Sort the values
    quicksort( X, I );
    // Delete duplicate entries
    size_t pos = 1;
    for ( size_t i = 1; i < X.size(); i++ ) {
        if ( X[i] != X[pos - 1] ) {
            X[pos] = X[i];
            I[pos] = I[i];
            pos++;
        }
    }
    X.resize( pos );
    I.resize( pos );
    // Fill the index array J
    for ( size_t i = 0; i < I.size(); i++ )
        J[I[i]] = i;
    for ( size_t i = 0; i < J.size(); i++ ) {
        if ( J[i] == neg_one )
            J[i] = J[i - 1];
    }
}


/************************************************************************
 * Subroutine to find the first element in X which is greater than Y     *
 * using a simple hashing technique.  This is the a faster method, but   *
 * requires the vector X to be in ascending order.                       *
 * Returns -1 if no value is larger.                                     *
 ************************************************************************/
template<class T>
size_t AMP::Utilities::findfirst( size_t n, const T *x, const T &value )
{
    AMP_INSIST( n > 0, "x must not be empty" );
    // Check if value is within the range of x
    if ( value <= x[0] )
        return 0;
    else if ( value > x[n - 1] )
        return n;
    // Perform the search
    size_t lower = 0;
    size_t upper = n - 1;
    size_t index;
    while ( ( upper - lower ) != 1 ) {
        index = ( upper + lower ) / 2;
        if ( x[index] >= value )
            upper = index;
        else
            lower = index;
    }
    index = upper;
    return index;
}

#endif
