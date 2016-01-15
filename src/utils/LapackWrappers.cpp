#ifndef included_LapackWrappers_tests_hpp_
#define included_LapackWrappers_tests_hpp_

#define NOMINMAX
#include "utils/LapackWrappers.h"
#include "utils/Utilities.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string.h>


// Choose the OS
#if defined( WIN32 ) || defined( _WIN32 ) || defined( WIN64 ) || defined( _WIN64 )
// Using windows
#define WINDOWS
#include <process.h>
#include <stdlib.h>
#include <windows.h>
#else
// Using some other operating system (probably linux)
#define LINUX
#include <pthread.h>
#include <unistd.h>
#endif


// Define some sizes of the problems
#define TEST_SIZE_VEC 10000  // Vector tests O(N)
#define TEST_SIZE_MATVEC 500 // Matrix-vector tests O(N^2)
#define TEST_SIZE_MAT 100    // Matrix-matrix / Dense solves tests O(N^3)
#define TEST_SIZE_TRI 1000   // Tridiagonal/banded tests


namespace AMP {


// Declare the individual tests
template <typename T>
static bool test_copy( int N, T &error );
template <typename T>
static bool test_scal( int N, T &error );
template <typename T>
static bool test_nrm2( int N, T &error );
template <typename T>
static bool test_asum( int N, T &error );
template <typename T>
static bool test_dot( int N, T &error );
template <typename T>
static bool test_axpy( int N, T &error );
template <typename T>
static bool test_gemv( int N, T &error );
template <typename T>
static bool test_gemm( int N, T &error );
template <typename T>
static bool test_gesv( int N, T &error );
template <typename T>
static bool test_gtsv( int N, T &error );
template <typename T>
static bool test_gbsv( int N, T &error );
template <typename T>
static bool test_getrf( int N, T &error );
template <typename T>
static bool test_gttrf( int N, T &error );
template <typename T>
static bool test_gbtrf( int N, T &error );
template <typename T>
static bool test_getrs( int N, T &error );
template <typename T>
static bool test_gttrs( int N, T &error );
template <typename T>
static bool test_gbtrs( int N, T &error );
template <typename T>
static bool test_getri( int N, T &error );


// Run a given test
template <>
int Lapack<double>::run_test( const char *routine, int N, double &error )
{
    std::string name( routine );
    // std::transform(name.begin(),name.end(),name.begin(),::tolower);
    int N_errors = 0;
    if ( name == "dcopy" ) {
        N_errors += test_copy<double>( N, error ) ? 1 : 0;
    } else if ( name == "dscal" ) {
        N_errors += test_scal<double>( N, error ) ? 1 : 0;
    } else if ( name == "dnrm2" ) {
        N_errors += test_nrm2<double>( N, error ) ? 1 : 0;
    } else if ( name == "daxpy" ) {
        N_errors += test_axpy<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgemv" ) {
        N_errors += test_gemv<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgemm" ) {
        N_errors += test_gemm<double>( N, error ) ? 1 : 0;
    } else if ( name == "dasum" ) {
        N_errors += test_asum<double>( N, error ) ? 1 : 0;
    } else if ( name == "ddot" ) {
        N_errors += test_dot<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgesv" ) {
        N_errors += test_gesv<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgtsv" ) {
        N_errors += test_gtsv<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgbsv" ) {
        N_errors += test_gbsv<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgetrf" ) {
        N_errors += test_getrf<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgttrf" ) {
        N_errors += test_gttrf<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgbtrf" ) {
        N_errors += test_gbtrf<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgetrs" ) {
        N_errors += test_getrs<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgttrs" ) {
        N_errors += test_gttrs<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgbtrs" ) {
        N_errors += test_gbtrs<double>( N, error ) ? 1 : 0;
    } else if ( name == "dgetri" ) {
        N_errors += test_getri<double>( N, error ) ? 1 : 0;
    } else {
        std::cerr << "Unknown test\n";
        return -1;
    }
    return N_errors;
}

// Run a given test
template <>
int Lapack<float>::run_test( const char *routine, int N, float &error )
{
    std::string name( routine );
    // std::transform(name.begin(),name.end(),name.begin(),::tolower);
    int N_errors = 0;
    if ( name == "scopy" ) {
        N_errors += test_copy<float>( N, error ) ? 1 : 0;
    } else if ( name == "sscal" ) {
        N_errors += test_scal<float>( N, error ) ? 1 : 0;
    } else if ( name == "snrm2" ) {
        N_errors += test_nrm2<float>( N, error ) ? 1 : 0;
    } else if ( name == "saxpy" ) {
        N_errors += test_axpy<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgemv" ) {
        N_errors += test_gemv<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgemm" ) {
        N_errors += test_gemm<float>( N, error ) ? 1 : 0;
    } else if ( name == "sasum" ) {
        N_errors += test_asum<float>( N, error ) ? 1 : 0;
    } else if ( name == "sdot" ) {
        N_errors += test_dot<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgesv" ) {
        N_errors += test_gesv<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgtsv" ) {
        N_errors += test_gtsv<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgbsv" ) {
        N_errors += test_gbsv<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgetrf" ) {
        N_errors += test_getrf<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgttrf" ) {
        N_errors += test_gttrf<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgbtrf" ) {
        N_errors += test_gbtrf<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgetrs" ) {
        N_errors += test_getrs<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgttrs" ) {
        N_errors += test_gttrs<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgbtrs" ) {
        N_errors += test_gbtrs<float>( N, error ) ? 1 : 0;
    } else if ( name == "sgetri" ) {
        N_errors += test_getri<float>( N, error ) ? 1 : 0;
    } else {
        std::cerr << "Unknown test\n";
        return -1;
    }
    return N_errors;
}

// Run all the tests
template <typename T>
int Lapack<T>::run_all_test()
{
    int N_errors = 0;
    int N        = 2; // We want two iterations to enure the test works for N>1
    T error;
    // Basic blas operations
    if ( test_copy<T>( N, error ) ) {
        printf( "test_copy failed (%e)\n", error );
        N_errors++;
    }
    if ( test_nrm2<T>( N, error ) ) {
        printf( "test_nrm2 failed (%e)\n", error );
        N_errors++;
    }
    if ( test_axpy<T>( N, error ) ) {
        printf( "test_axpy failed (%e)\n", error );
        N_errors++;
    }
    if ( test_asum<T>( N, error ) ) {
        printf( "test_asum failed (%e)\n", error );
        N_errors++;
    }
    if ( test_dot<T>( N, error ) ) {
        printf( "test_dot failed (%e)\n", error );
        N_errors++;
    }
    // Matrix blas operations
    if ( test_gemv<T>( N, error ) ) {
        printf( "test_gemv failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gemm<T>( N, error ) ) {
        printf( "test_gemm failed (%e)\n", error );
        N_errors++;
    }
    // Linear solves
    if ( test_gesv<T>( N, error ) ) {
        printf( "test_gesv failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gtsv<T>( N, error ) ) {
        printf( "test_gtsv failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gbsv<T>( N, error ) ) {
        printf( "test_gbsv failed (%e)\n", error );
        N_errors++;
    }
    // Linear factorizations
    if ( test_getrf<T>( N, error ) ) {
        printf( "test_getrf failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gttrf<T>( N, error ) ) {
        printf( "test_gttrf failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gbtrf<T>( N, error ) ) {
        printf( "test_gbtrf failed (%e)\n", error );
        N_errors++;
    }
    // Solve using factorization
    if ( test_getrs<T>( N, error ) ) {
        printf( "test_getrs failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gttrs<T>( N, error ) ) {
        printf( "test_gttrs failed (%e)\n", error );
        N_errors++;
    }
    if ( test_gbtrs<T>( N, error ) ) {
        printf( "test_gbtrs failed (%e)\n", error );
        N_errors++;
    }
    // Inverse using factorization
    if ( test_getri<T>( N, error ) ) {
        printf( "test_getri failed (%e)\n", error );
        N_errors++;
    }
    return N_errors > 0;
}


// Fill a vector with random T precision data
template <typename T>
static inline void random( int N, T *data )
{
    const T rmax = static_cast<T>( RAND_MAX - 1 );
    for ( int i = 0; i < N; i++ )
        data[i] = static_cast<T>( rand() / rmax + 1e-7 * static_cast<T>( rand() ) / rmax );
}


// Check if two vectors are approximately equal
template <typename T>
static inline bool approx_equal( int N, const T *x1, const T *x2, const T tol = 1e-12 )
{
    bool pass = true;
    for ( int i = 0; i < N; i++ )
        pass = pass && fabs( x1[i] - x2[i] ) <= tol * 0.5 * fabs( x1[i] + x2[i] );
    return pass;
}


// Return the L2 norm
template <typename T>
static inline T L2Norm( int N, const T *x )
{
    T norm = 0.0;
    for ( int i = 0; i < N; i++ )
        norm += x[i] * x[i];
    return sqrt( norm );
}


// Return the L2 error
template <typename T>
static inline T L2Error( int N, const T *x1, const T *x2 )
{
    T norm = 0.0;
    for ( int i = 0; i < N; i++ )
        norm += ( x1[i] - x2[i] ) * ( x1[i] - x2[i] );
    return sqrt( norm );
}


// Test dcopy
template <typename T>
static bool test_copy( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    T *x1       = new T[K];
    T *x2       = new T[K];
    random( K, x1 );
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        memset( x2, 0xB6, K * sizeof( T ) );
        Lapack<T>::copy( K, x1, 1, x2, 1 );
        if ( !approx_equal( K, x1, x2 ) )
            N_errors++;
    }
    delete[] x1;
    delete[] x2;
    return N_errors > 0;
}

// Test dcopy
template <typename T>
static bool test_scal( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    T *x0       = new T[K];
    T *x1       = new T[K];
    T *x2       = new T[K];
    random( K, x0 );
    const T pi = static_cast<T>( 3.141592653589793 );
    for ( int j  = 0; j < K; j++ )
        x1[j]    = pi * x0[j];
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, x1, K * sizeof( T ) );
        Lapack<T>::scal( K, pi, x0, 1 );
        if ( !approx_equal( K, x1, x2, 10 * std::numeric_limits<T>::epsilon() ) )
            N_errors++;
    }
    delete[] x0;
    delete[] x1;
    delete[] x2;
    return N_errors > 0;
}

// Test dnrm2
template <typename T>
static bool test_nrm2( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    T *x        = new T[K];
    random( K, x );
    T ans1 = 0.0;
    for ( int j = 0; j < K; j++ )
        ans1 += x[j] * x[j];
    ans1         = sqrt( ans1 );
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        T ans2 = Lapack<T>::nrm2( K, x, 1 );
        error  = std::max( error, std::abs( ans1 - ans2 ) );
        if ( std::abs( ans1 - ans2 ) > K * std::numeric_limits<T>::epsilon() )
            N_errors++;
    }
    delete[] x;
    return N_errors > 0;
}

// Test dasum
template <typename T>
static bool test_asum( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    // Maximum roundoff error that is acceptible is determined by the
    //    error from adding a series of numbers from (0,1).
    //    The expected error is much less than the maximum error.
    const T max_error = K * K / 2 * std::numeric_limits<T>::epsilon();
    // Create a random set of numbers and the sum (simple method)
    T *x = new T[K];
    random( K, x );
    T ans1 = 0;
    for ( int j = 0; j < K; j++ )
        ans1 += std::abs( x[j] );
    // Check dasum
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        T ans2 = Lapack<T>::asum( K, x, 1 );
        error  = std::max( error, std::abs( ( ans1 - ans2 ) / K ) );
        if ( std::abs( ans1 - ans2 ) > max_error )
            N_errors++;
    }
    delete[] x;
    return N_errors > 0;
}

// Test ddot
template <typename T>
static bool test_dot( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    T *x1       = new T[K];
    T *x2       = new T[K];
    random( K, x1 );
    random( K, x2 );
    T ans1 = 0.0;
    for ( int j = 0; j < K; j++ )
        ans1 += x1[j] * x2[j];
    int N_errors = 0;
    error        = 0;
    for ( int i = 0; i < N; i++ ) {
        T ans2 = Lapack<T>::dot( K, x1, 1, x2, 1 );
        error  = std::max( error, std::abs( ans1 - ans2 ) / K );
        if ( std::abs( ans1 - ans2 ) > 10 * K * std::numeric_limits<T>::epsilon() )
            N_errors++;
    }
    delete[] x1;
    delete[] x2;
    return N_errors > 0;
}

// Test daxpy
template <typename T>
static bool test_axpy( int N, T &error )
{
    const int K = TEST_SIZE_VEC;
    T *x        = new T[K];
    T *y0       = new T[K];
    T *y1       = new T[K];
    T *y2       = new T[K];
    random( K, x );
    random( K, y0 );
    const T pi = static_cast<T>( 3.141592653589793 );
    for ( int j = 0; j < K; j++ )
        y1[j]   = y0[j] + pi * x[j];
    error       = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( y2, y0, K * sizeof( T ) );
        Lapack<T>::axpy( K, pi, x, 1, y2, 1 );
        T err = L2Error( K, y1, y2 );
        error = std::max( error, err );
    }
    bool fail = error > std::numeric_limits<T>::epsilon();
    NULL_USE( y1 );
    delete[] x;
    delete[] y0;
    delete[] y1;
    delete[] y2;
    return fail;
}

// Test dgemv
template <typename T>
static bool test_gemv( int N, T &error )
{
    const int K = TEST_SIZE_MATVEC;
    T *A        = new T[K * K];
    T *x        = new T[K];
    T *y        = new T[K];
    T *y1       = new T[K];
    T *y2       = new T[K];
    random( K * K, A );
    random( K, x );
    random( K, y );
    const T alpha = static_cast<T>( 3.141592653589793 );
    const T beta  = static_cast<T>( 1.414213562373095 );
    for ( int j = 0; j < K; j++ ) {
        y1[j] = beta * y[j];
        for ( int k = 0; k < K; k++ )
            y1[j] += alpha * A[j + k * K] * x[k];
    }
    int N_errors = 0;
    error        = 0;
    T norm       = L2Norm( K, y1 );
    for ( int i = 0; i < N; i++ ) {
        memcpy( y2, y, K * sizeof( T ) );
        Lapack<T>::gemv( 'N', K, K, alpha, A, K, x, 1, beta, y2, 1 );
        error = std::max( error, L2Error( K, y1, y2 ) / norm );
        if ( !approx_equal( K, y1, y2, K * std::numeric_limits<T>::epsilon() ) )
            N_errors++;
    }
    delete[] A;
    delete[] x;
    delete[] y;
    delete[] y1;
    delete[] y2;
    return N_errors > 0;
}

// Test dgemm
template <typename T>
static bool test_gemm( int N, T &error )
{
    const int K = TEST_SIZE_MAT;
    T *A        = new T[K * K];
    T *B        = new T[K * K];
    T *C        = new T[K * K];
    T *C1       = new T[K * K];
    T *C2       = new T[K * K];
    random( K * K, A );
    random( K * K, B );
    random( K * K, C );
    const T alpha = static_cast<T>( 3.141592653589793 );
    const T beta  = static_cast<T>( 1.414213562373095 );
    for ( int i = 0; i < K; i++ ) {
        for ( int j = 0; j < K; j++ ) {
            C1[i + j * K] = beta * C[i + j * K];
            for ( int k = 0; k < K; k++ )
                C1[i + j * K] += alpha * A[i + k * K] * B[k + j * K];
        }
    }
    int N_errors = 0;
    error        = 0;
    T norm       = L2Norm( K * K, C1 );
    for ( int i = 0; i < N; i++ ) {
        memcpy( C2, C, K * K * sizeof( T ) );
        Lapack<T>::gemm( 'N', 'N', K, K, K, alpha, A, K, B, K, beta, C2, K );
        error = std::max( error, L2Error( K * K, C1, C2 ) / norm );
        if ( !approx_equal( K * K, C1, C2, K * 10 * std::numeric_limits<T>::epsilon() ) )
            N_errors++;
    }
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] C1;
    delete[] C2;
    return N_errors > 0;
}

// Test dgesv
template <typename T>
static bool test_gesv( int N, T &error )
{
    // Test solving a diagonal matrix
    const int K = TEST_SIZE_MAT;
    T *A        = new T[K * K];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    memset( A, 0, K * K * sizeof( T ) );
    random( K, x2 );
    random( K, b );
    for ( int k = 0; k < K; k++ ) {
        A[k + k * K] = x2[k] + std::numeric_limits<T>::epsilon();
        x1[k]        = b[k] / ( x2[k] + std::numeric_limits<T>::epsilon() );
    }
    int N_errors = 0;
    error        = 0;
    T norm       = L2Norm( K, x1 );
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( T ) );
        int err = 0;
        Lapack<T>::gesv( K, 1, A, K, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        error = std::max( error, L2Error( K, x1, x2 ) / norm );
        if ( !approx_equal( K, x1, x2, K * 10 * std::numeric_limits<T>::epsilon() ) )
            N_errors++;
    }
    delete[] A;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgtsv
template <typename T>
static bool test_gtsv( int N, T &error )
{
    // Test solving a tri-diagonal matrix by comparing to dgtsv
    const int K = TEST_SIZE_TRI;
    T *A        = new T[K * K];
    T *D        = new T[K];
    T *D2       = new T[K];
    T *DL       = new T[K - 1];
    T *DL2      = new T[K - 1];
    T *DU       = new T[K - 1];
    T *DU2      = new T[K - 1];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    memset( A, 0, K * K * sizeof( T ) );
    random( K, D );
    random( K - 1, DL );
    random( K - 1, DU );
    random( K, b );
    A[0] = D[0];
    for ( int k = 1; k < K; k++ ) {
        A[k + k * K]         = D[k];
        A[( k - 1 ) + k * K] = DU[k - 1];
        A[k + ( k - 1 ) * K] = DL[k - 1];
    }
    memcpy( x1, b, K * sizeof( T ) );
    int err = 0;
    Lapack<T>::gesv( K, 1, A, K, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( T ) );
        memcpy( D2, D, K * sizeof( T ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( T ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( T ) );
        Lapack<T>::gtsv( K, 1, DL2, D2, DU2, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        if ( err != 0 )
            printf( "Error calling gtsv (%i)\n", err );
        T err2 = L2Error( N, x1, x2 );
        T norm = L2Norm( N, x1 );
        error  = std::max( error, err2 / norm );
    }
    const T tol = static_cast<T>( 2e4 ) * std::numeric_limits<T>::epsilon();
    if ( error > tol ) {
        printf( "test_gtsv error (%e) exceeded tolerance (%e)\n", error, tol );
        N_errors++;
    }
    delete[] A;
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}
// Test dgtsv
template <typename T>
static bool test_gbsv( int N, T &error )
{
    // Test solving a banded-diagonal matrix by comparing to dgtsv
    //    N = 6, KL = 2, KU = 1:
    //        *    *    *    +    +    +
    //        *    *    +    +    +    +
    //        *   a12  a23  a34  a45  a56
    //       a11  a22  a33  a44  a55  a66
    //       a21  a32  a43  a54  a65   *
    //       a31  a42  a53  a64   *    *
    const int K  = TEST_SIZE_TRI;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    T *A         = new T[K * K];
    T *AB        = new T[K * K2];
    T *AB2       = new T[K * K2];
    T *x1        = new T[K];
    T *x2        = new T[K];
    T *b         = new T[K];
    int *IPIV    = new int[K];
    random( K * K2, AB );
    random( K, b );
    memset( A, 0, K * K * sizeof( T ) );
    for ( int k = 0; k < K; k++ ) {
        for ( int k2 = -KL; k2 <= KU; k2++ ) {
            if ( k + k2 < 0 || k + k2 >= K ) {
                continue;
            }
            A[k + k2 + k * K] = AB[k2 + 2 * KL + k * K2];
        }
    }
    memcpy( x1, b, K * sizeof( T ) );
    int err = 0;
    Lapack<T>::gesv( K, 1, A, K, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( x2, b, K * sizeof( T ) );
        memcpy( AB2, AB, K * K2 * sizeof( T ) );
        Lapack<T>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        T norm = L2Norm( K, x1 );
        T err2 = L2Error( K, x1, x2 );
        error  = std::max( error, err2 / norm );
    }
    const double tol = 2000.0 * std::numeric_limits<T>::epsilon();
    if ( error > tol ) {
        printf( "test_gbsv error (%e) exceeded tolerance (%e)\n", error, tol );
        N_errors++;
    }
    delete[] A;
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgetrf
template <typename T>
static bool test_getrf( int N, T &error )
{
    // Check dgetrf by performing a factorization and solve and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    T *A        = new T[K * K];
    T *A2       = new T[K * K];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    random( K * K, A );
    random( K, b );
    memcpy( A2, A, K * K * sizeof( T ) );
    memcpy( x1, b, K * sizeof( T ) );
    int err = 0;
    Lapack<T>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( A2, A, K * K * sizeof( T ) );
        Lapack<T>::getrf( K, K, A2, K, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( T ) );
    Lapack<T>::getrs( 'N', K, 1, A2, K, IPIV, x2, K, err );
    T norm = L2Norm( K, x1 );
    T err2 = L2Error( K, x1, x2 );
    if ( err2 > 10.0 * norm * std::numeric_limits<T>::epsilon() )
        N_errors++;
    error = err2 / norm;
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgttrf
template <typename T>
static bool test_gttrf( int N, T &error )
{
    // Check dgttrf by performing a factorization and solve and comparing to dgtsv
    const int K = TEST_SIZE_TRI;
    T *D        = new T[K];
    T *D2       = new T[K];
    T *DL       = new T[K - 1];
    T *DL2      = new T[K - 1];
    T *DU       = new T[K - 1];
    T *DU2      = new T[K - 1];
    T *DU3      = new T[K - 2];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    random( K, D );
    random( K - 1, DL );
    random( K - 1, DU );
    random( K, b );
    memcpy( x1, b, K * sizeof( T ) );
    memcpy( D2, D, K * sizeof( T ) );
    memcpy( DL2, DL, ( K - 1 ) * sizeof( T ) );
    memcpy( DU2, DU, ( K - 1 ) * sizeof( T ) );
    int err = 0;
    Lapack<T>::gtsv( K, 1, DL2, D2, DU2, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( D2, D, K * sizeof( T ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( T ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( T ) );
        Lapack<T>::gttrf( K, DL2, D2, DU2, DU3, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( T ) );
    Lapack<T>::gttrs( 'N', K, 1, DL2, D2, DU2, DU3, IPIV, x2, K, err );
    T norm = L2Norm( K, x1 );
    T err2 = L2Error( K, x1, x2 );
    if ( err2 > 10.0 * norm * std::numeric_limits<T>::epsilon() )
        N_errors++;
    error = err2 / norm;
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] DU3;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgbtrf
template <typename T>
static bool test_gbtrf( int N, T &error )
{
    // Check dgbtrf by performing a factorization and solve and comparing to dgbsv
    const int K  = TEST_SIZE_TRI;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    T *AB        = new T[K * K2];
    T *AB2       = new T[K * K2];
    T *x1        = new T[K];
    T *x2        = new T[K];
    T *b         = new T[K];
    int *IPIV    = new int[K];
    random( K * K2, AB );
    random( K, b );
    int err = 0;
    memcpy( x1, b, K * sizeof( T ) );
    memcpy( AB2, AB, K * K2 * sizeof( T ) );
    Lapack<T>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x1, K, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( AB2, AB, K * K2 * sizeof( T ) );
        Lapack<T>::gbtrf( K, K, KL, KU, AB2, K2, IPIV, err );
        N_errors += err == 0 ? 0 : 1;
    }
    memcpy( x2, b, K * sizeof( T ) );
    Lapack<T>::gbtrs( 'N', K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
    T norm = L2Norm( K, x1 );
    T err2 = L2Error( K, x1, x2 );
    if ( err2 > 10.0 * norm * std::numeric_limits<T>::epsilon() )
        N_errors++;
    error = err2 / norm;
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgetrs
template <typename T>
static bool test_getrs( int N, T &error )
{
    // Check dgetrs by performing a factorization and solve and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    T *A        = new T[K * K];
    T *A2       = new T[K * K];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    random( K * K, A );
    random( K, b );
    int err = 0;
    memcpy( A2, A, K * K * sizeof( T ) );
    memcpy( x1, b, K * sizeof( T ) );
    Lapack<T>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    int N_errors = 0;
    Lapack<T>::getrf( K, K, A, K, IPIV, err );
    for ( int i = 0; i < N; i++ ) {
        memcpy( A2, A, K * K * sizeof( T ) );
        memcpy( x2, b, K * sizeof( T ) );
        Lapack<T>::getrs( 'N', K, 1, A2, K, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        T norm = L2Norm( K, x1 );
        T err2 = L2Error( K, x1, x2 );
        if ( err > 10.0 * norm * std::numeric_limits<T>::epsilon() )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgttrs
template <typename T>
static bool test_gttrs( int N, T &error )
{
    // Check dgttrs by performing a factorization and solve and comparing to dgtsv
    const int K = TEST_SIZE_TRI;
    T *D        = new T[K];
    T *D2       = new T[K];
    T *DL       = new T[K - 1];
    T *DL2      = new T[K - 1];
    T *DU       = new T[K - 1];
    T *DU2      = new T[K - 1];
    T *DU3      = new T[K - 2];
    T *DU4      = new T[K - 2];
    T *x1       = new T[K];
    T *x2       = new T[K];
    T *b        = new T[K];
    int *IPIV   = new int[K];
    random( K, D );
    random( K - 1, DL );
    random( K - 1, DU );
    random( K, b );
    int err = 0;
    memcpy( x1, b, K * sizeof( T ) );
    memcpy( D2, D, K * sizeof( T ) );
    memcpy( DL2, DL, ( K - 1 ) * sizeof( T ) );
    memcpy( DU2, DU, ( K - 1 ) * sizeof( T ) );
    Lapack<T>::gtsv( K, 1, DL2, D2, DU2, x1, K, err );
    Lapack<T>::gttrf( K, DL, D, DU, DU3, IPIV, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( D2, D, K * sizeof( T ) );
        memcpy( DL2, DL, ( K - 1 ) * sizeof( T ) );
        memcpy( DU2, DU, ( K - 1 ) * sizeof( T ) );
        memcpy( DU4, DU3, ( K - 2 ) * sizeof( T ) );
        memcpy( x2, b, K * sizeof( T ) );
        Lapack<T>::gttrs( 'N', K, 1, DL2, D2, DU2, DU4, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        T norm = L2Norm( K, x1 );
        T err2 = L2Error( K, x1, x2 );
        if ( err2 > 10.0 * norm * std::numeric_limits<T>::epsilon() )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] D;
    delete[] D2;
    delete[] DL;
    delete[] DL2;
    delete[] DU;
    delete[] DU2;
    delete[] DU3;
    delete[] DU4;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgbtrs
template <typename T>
static bool test_gbtrs( int N, T &error )
{
    // Check dgbtrs by performing a factorization and solve and comparing to dgbsv
    const int K  = TEST_SIZE_TRI;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2 * KL + KU + 1;
    T *AB        = new T[K * K2];
    T *AB2       = new T[K * K2];
    T *x1        = new T[K];
    T *x2        = new T[K];
    T *b         = new T[K];
    int *IPIV    = new int[K];
    random( K * K2, AB );
    random( K, b );
    int err = 0;
    memcpy( x1, b, K * sizeof( T ) );
    memcpy( AB2, AB, K * K2 * sizeof( T ) );
    Lapack<T>::gbsv( K, KL, KU, 1, AB2, K2, IPIV, x1, K, err );
    Lapack<T>::gbtrf( K, K, KL, KU, AB, K2, IPIV, err );
    int N_errors = 0;
    for ( int i = 0; i < N; i++ ) {
        memcpy( AB2, AB, K * K2 * sizeof( T ) );
        memcpy( x2, b, K * sizeof( T ) );
        Lapack<T>::gbtrs( 'N', K, KL, KU, 1, AB2, K2, IPIV, x2, K, err );
        N_errors += err == 0 ? 0 : 1;
        T norm = L2Norm( K, x1 );
        T err2 = L2Error( K, x1, x2 );
        if ( err2 > 10.0 * norm * std::numeric_limits<T>::epsilon() )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] AB;
    delete[] AB2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    return N_errors > 0;
}

// Test dgetri
template <typename T>
static bool test_getri( int N, T &error )
{
    // Check dgetri by performing a factorization, calculating the inverse,
    //   multiplying the rhs, and comparing to dgesv
    const int K     = TEST_SIZE_MAT;
    const int LWORK = 8 * K;
    T *A            = new T[K * K];
    T *A2           = new T[K * K];
    T *x1           = new T[K];
    T *x2           = new T[K];
    T *b            = new T[K];
    int *IPIV       = new int[K];
    T *WORK         = new T[LWORK];
    random( K * K, A );
    random( K, b );
    int err = 0;
    memcpy( A2, A, K * K * sizeof( T ) );
    memcpy( x1, b, K * sizeof( T ) );
    Lapack<T>::gesv( K, 1, A2, K, IPIV, x1, K, err );
    int N_errors = 0;
    Lapack<T>::getrf( K, K, A, K, IPIV, err );
    for ( int i = 0; i < N; i++ ) {
        // Compute the inverse
        memcpy( A2, A, K * K * sizeof( T ) );
        Lapack<T>::getri( K, A2, K, IPIV, WORK, LWORK, err );
        N_errors += err == 0 ? 0 : 1;
        // Perform the mat-vec
        memset( x2, 0xB6, K * sizeof( T ) );
        Lapack<T>::gemv( 'N', K, K, 1, A2, K, b, 1, 0, x2, 1 );
        // Check the result
        T norm = L2Norm( K, x1 );
        T err2 = L2Error( K, x1, x2 );
        if ( err2 > 100 * norm * std::numeric_limits<T>::epsilon() )
            N_errors++;
        error = std::max( error, err2 / norm );
    }
    delete[] A;
    delete[] A2;
    delete[] x1;
    delete[] x2;
    delete[] b;
    delete[] IPIV;
    delete[] WORK;
    return N_errors > 0;
}


/******************************************************************
* Print all of the machine parameters by lamch                   *
******************************************************************/
template <typename T>
void Lapack<T>::print_machine_parameters()
{
    printf( "eps   = %13.6e    relative machine precision\n", Lapack<T>::lamch( 'E' ) );
    printf( "sfmin = %13.6e    safe minimum\n", Lapack<T>::lamch( 'S' ) );
    printf( "base  = %13.6e    base of the machine\n", Lapack<T>::lamch( 'B' ) );
    printf( "prec  = %13.6e    eps*base\n", Lapack<T>::lamch( 'P' ) );
    printf( "t     = %13.6e    number of digits in the mantissa\n", Lapack<T>::lamch( 'N' ) );
    printf( "rnd   = %13.6e    1.0 when rounding occurs in addition, 0.0 otherwise\n",
            Lapack<T>::lamch( 'R' ) );
    printf( "emin  = %13.6e    minimum exponent before underflow\n", Lapack<T>::lamch( 'M' ) );
    printf( "rmin  = %13.6e    underflow threshold - base**(emin-1)\n", Lapack<T>::lamch( 'U' ) );
    printf( "emax  = %13.6e    largest exponent before overflow\n", Lapack<T>::lamch( 'L' ) );
    printf( "rmax  = %13.6e    overflow threshold - (base**emax)*(1-eps)\n",
            Lapack<T>::lamch( 'O' ) );
}


/******************************************************************
* Some inline functions to acquire/release a mutex                *
******************************************************************/
#ifdef WINDOWS
static HANDLE LapackWrappers_lock_queue = CreateMutex( NULL, FALSE, NULL );
#elif defined( LINUX )
static pthread_mutex_t LapackWrappers_lock_queue;
static int LapackWrappers_lock_queue_error = pthread_mutex_init( &LapackWrappers_lock_queue, NULL );
#else
#error Not programmed
#endif
#ifdef WINDOWS
template <typename T>
void Lapack<T>::get_lock()
{
    WaitForSingleObject( LapackWrappers_lock_queue, INFINITE );
}
#elif defined( LINUX )
template <typename T>
void Lapack<T>::get_lock()
{
    int error = pthread_mutex_lock( &LapackWrappers_lock_queue );
    if ( error == -1 )
        printf( "Error locking mutex" );
}
#else
#error Not programmed
#endif
#ifdef WINDOWS
template <typename T>
void Lapack<T>::release_lock()
{
    bool success = ReleaseMutex( LapackWrappers_lock_queue ) != 0;
    if ( !success )
        printf( "Error unlocking mutex" );
}
#elif defined( LINUX )
template <typename T>
void Lapack<T>::release_lock()
{
    int error = pthread_mutex_unlock( &LapackWrappers_lock_queue );
    if ( error == -1 )
        printf( "Error unlocking mutex" );
}
#else
#error Not programmed
#endif


} // namespace


// Explicit instantiations
template class AMP::Lapack<double>;
template class AMP::Lapack<float>;


#endif
