#include "utils/LapackWrappers.h"
#include <iostream>
#include <cstdio>
#include <limits>
#include <string.h>
#include <algorithm>
#include <cmath>


using namespace AMP;


// Choose the OS 
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    // Using windows
    #define WINDOWS
    #include <stdlib.h>
    #include <windows.h>
    #include <process.h>
#else
    // Using some other operating system (probably linux)
    #define LINUX
    #include <pthread.h>
    #include <unistd.h>
#endif


// Declare the individual tests
static bool test_dcopy( int N );
static bool test_dscal( int N );
static bool test_dnrm2( int N );
static bool test_dasum( int N );
static bool test_ddot(  int N );
static bool test_daxpy( int N );
static bool test_dgemv( int N );
static bool test_dgemm( int N );
static bool test_dgesv( int N );
static bool test_dgtsv( int N );
static bool test_dgbsv( int N );
static bool test_dgetrf( int N );
static bool test_dgttrf( int N );
static bool test_dgbtrf( int N );
static bool test_dgetrs( int N );
static bool test_dgttrs( int N );
static bool test_dgbtrs( int N );
static bool test_dgetri( int N );


// Define some sizes of the problems
#define TEST_SIZE_VEC  10000    // Vector tests
#define TEST_SIZE_MAT  100      // Matrix tests
#define TEST_SIZE_tri  500      // Tridiagonal/banded tests

// Run a given test
int AMP::Lapack::run_test( const char* routine, int N )
{
    std::string name(routine);
    //std::transform(name.begin(),name.end(),name.begin(),::tolower);
    int N_errors = 0;
    if ( name == "dcopy" ) {
        N_errors += test_dcopy( N ) ? 1:0;
    } else if ( name == "dscal" ) {
        N_errors += test_dscal(N) ? 1 : 0;
    } else if ( name == "dnrm2" ) {
        N_errors += test_dnrm2(N) ? 1 : 0;
    } else if ( name == "daxpy" ) {
        N_errors += test_daxpy(N) ? 1 : 0;
    } else if ( name == "dgemv" ) {
        N_errors += test_dgemv(N) ? 1 : 0;
    } else if ( name == "dgemm" ) {
        N_errors += test_dgemm(N) ? 1 : 0;
    } else if ( name == "dasum" ) {
        N_errors += test_dasum(N) ? 1 : 0;
    } else if ( name == "ddot" ) {
        N_errors += test_ddot(N)  ? 1 : 0;
    } else if ( name == "dgesv" ) {
        N_errors += test_dgesv(N) ? 1 : 0;
    } else if ( name == "dgtsv" ) {
        N_errors += test_dgtsv(N) ? 1 : 0;
    } else if ( name == "dgbsv" ) {
        N_errors += test_dgbsv(N) ? 1 : 0;
    } else if ( name == "dgetrf" ) {
        N_errors += test_dgetrf(N) ? 1 : 0;
    } else if ( name == "dgttrf" ) {
        N_errors += test_dgttrf(N) ? 1 : 0;
    } else if ( name == "dgbtrf" ) {
        N_errors += test_dgbtrf(N) ? 1 : 0;
    } else if ( name == "dgetrs" ) {
        N_errors += test_dgetrs(N) ? 1 : 0;
    } else if ( name == "dgttrs" ) {
        N_errors += test_dgttrs(N) ? 1 : 0;
    } else if ( name == "dgbtrs" ) {
        N_errors += test_dgbtrs(N) ? 1 : 0;
    } else if ( name == "dgetri" ) {
        N_errors += test_dgetri(N) ? 1 : 0;
    } else { 
        std::cerr << "Unknown test\n";
        return -1;
    }
    return N_errors;
}


// Run all the tests
int AMP::Lapack::run_all_test( )
{
    int N_errors = 0;
    int N = 2;  // We want two iterations to enure the test works for N>1
    // Basic blas operations
    if ( test_dcopy(N)!=0 ) { printf("test_dcopy failed\n"); N_errors++; }
    if ( test_dnrm2(N)!=0 ) { printf("test_dnrm2 failed\n"); N_errors++; }
    if ( test_daxpy(N)!=0 ) { printf("test_daxpy failed\n"); N_errors++; }
    if ( test_dasum(N)!=0 ) { printf("test_dasum failed\n"); N_errors++; }
    if ( test_ddot(N)!=0  ) { printf("test_ddot failed\n");  N_errors++; }
    // Matrix blas operations
    if ( test_dgemv(N)!=0 ) { printf("test_dgemv failed\n"); N_errors++; }
    if ( test_dgemm(N)!=0 ) { printf("test_dgemm failed\n"); N_errors++; }
    // Linear solves
    if ( test_dgesv(N)!=0 ) { printf("test_dgesv failed\n"); N_errors++; }
    if ( test_dgtsv(N)!=0 ) { printf("test_dgtsv failed\n"); N_errors++; }
    if ( test_dgbsv(N)!=0 ) { printf("test_dgbsv failed\n"); N_errors++; }
    // Linear factorizations
    if ( test_dgetrf(N)!=0 ) { printf("test_dgetrf failed\n"); N_errors++; }
    if ( test_dgttrf(N)!=0 ) { printf("test_dgttrf failed\n"); N_errors++; }
    if ( test_dgbtrf(N)!=0 ) { printf("test_dgbtrf failed\n"); N_errors++; }
    // Solve using factorization
    if ( test_dgetrs(N)!=0 ) { printf("test_dgetrs failed\n"); N_errors++; }
    if ( test_dgttrs(N)!=0 ) { printf("test_dgttrs failed\n"); N_errors++; }
    if ( test_dgbtrs(N)!=0 ) { printf("test_dgbtrs failed\n"); N_errors++; }
    // Inverse using factorization
    if ( test_dgetri(N)!=0 ) { printf("test_dgetri failed\n"); N_errors++; }
    return N_errors>0;
}


// Fill a vector with random double precision data
static inline void random( int N, double *data )
{
    const double rmax = static_cast<double>(RAND_MAX-1);
    for (int i=0; i<N; i++)
        data[i] = static_cast<double>(rand())/rmax + 1e-7*static_cast<double>(rand())/rmax;
}


// Check if two vectors are approximately equal
static inline bool approx_equal( int N, const double *x1, const double *x2, const double tol = 1e-12 )
{
    bool pass = true;
    for (int i=0; i<N; i++)
        pass = pass && fabs(x1[i] - x2[i]) <= tol*0.5*fabs(x1[i] + x2[i]);
    return pass;
}


// Return the L2 norm
static inline double L2Norm( int N, const double *x )
{
    double norm = 0.0;
    for (int i=0; i<N; i++)
        norm += x[i]*x[i];
    return sqrt(norm);
}


// Return the L2 error
static inline double L2Error( int N, const double *x1, const double *x2 )
{
    double norm = 0.0;
    for (int i=0; i<N; i++)
        norm += (x1[i]-x2[i])*(x1[i]-x2[i]);
    return sqrt(norm);
}



// Test dcopy
static bool test_dcopy( int N )
{
    const int K = TEST_SIZE_VEC;
    double *x1 = new double[K];
    double *x2 = new double[K];
    random(K,x1);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memset(x2,0xB6,K*sizeof(double));
        Lapack::dcopy(K,x1,1,x2,1);
        if ( !approx_equal(K,x1,x2,0) )
            N_errors++;
    }
    delete [] x1;
    delete [] x2;
    return N_errors>0;
}

// Test dcopy
static bool test_dscal( int N )
{
    const int K = TEST_SIZE_VEC;
    double *x0 = new double[K];
    double *x1 = new double[K];
    double *x2 = new double[K];
    random(K,x0);
    const double pi = 3.141592653589793;
    for (int j=0; j<K; j++)
        x1[j] = pi*x0[j];
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(x2,x1,K*sizeof(double));
        Lapack::dscal(K,pi,x0,1);
        if ( !approx_equal(K,x1,x2,1e-14) )
            N_errors++;
    }
    delete [] x0;
    delete [] x1;
    delete [] x2;
    return N_errors>0;
}

// Test dnrm2
static bool test_dnrm2( int N )
{
    const int K = TEST_SIZE_VEC;
    double *x = new double[K];
    random(K,x);
    double ans1 = 0.0;
    for (int j=0; j<K; j++)
        ans1 += x[j]*x[j];
    ans1 = sqrt(ans1);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        double ans2 = Lapack::dnrm2(K,x,1);
        if ( fabs(ans1-ans2)>K*1e-15 )
            N_errors++;
    }
    delete [] x;
    return N_errors>0;
}

// Test dasum
static bool test_dasum( int N )
{
    const int K = TEST_SIZE_VEC;
    // Maximum roundoff error that is acceptible is determined by the
    //    error from adding a series of numbers from (0,1).
    //    The expected error is much less than the maximum error.
    const double max_error = K*K/2*std::numeric_limits<double>::epsilon();
    // Create a random set of numbers and the sum (simple method)
    double *x = new double[K];
    random(K,x);
    double ans1 = 0;
    for (int j=0; j<K; j++)
        ans1 += fabs(x[j]);
    // Check dasum
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        double ans2 = Lapack::dasum(K,x,1);
        if ( fabs(ans1-ans2)>max_error )
            N_errors++;
    }
    delete [] x;
    return N_errors>0;
}

// Test ddot
static bool test_ddot( int N )
{
    const int K = TEST_SIZE_VEC;
    double *x1 = new double[K];
    double *x2 = new double[K];
    random(K,x1);
    random(K,x2);
    double ans1 = 0.0;
    for (int j=0; j<K; j++)
        ans1 += x1[j]*x2[j];
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        double ans2 = Lapack::ddot(K,x1,1,x2,1);
        if ( fabs(ans1-ans2)>2*K*1e-15 )
            N_errors++;
    }
    delete [] x1;
    delete [] x2;
    return N_errors>0;
}

// Test daxpy
static bool test_daxpy( int N )
{
    const int K = TEST_SIZE_VEC;
    double *x = new double[K];
    double *y0 = new double[K];
    double *y1 = new double[K];
    double *y2 = new double[K];
    random(K,x);
    random(K,y1);
    const double pi = 3.141592653589793;
    for (int j=0; j<K; j++)
        y1[j] = y0[j] + pi*x[j];
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(y2,y0,K*sizeof(double));
        Lapack::daxpy(K,pi,x,1,y2,1);
        if ( !approx_equal(K,y1,y2,1e-14) )
            N_errors++;
    }
    delete [] x;
    delete [] y0;
    delete [] y1;
    delete [] y2;
    return N_errors>0;
}

// Test dgemv
static bool test_dgemv( int N )
{
    const int K = 500;
    double *A = new double[K*K];
    double *x = new double[K];
    double *y = new double[K];
    double *y1 = new double[K];
    double *y2 = new double[K];
    random(K*K,A);
    random(K,x);
    random(K,y);
    const double alpha = 3.141592653589793;
    const double beta  = 1.414213562373095;
    int N_errors = 0;
    for (int j=0; j<K; j++) {
        y1[j] = beta*y[j];
        for (int k=0; k<K; k++)
            y1[j] += alpha*A[j+k*K]*x[k];
    }
    for (int i=0; i<N; i++) {
        memcpy(y2,y,K*sizeof(double));
        Lapack::dgemv('N',K,K,alpha,A,K,x,1,beta,y2,1);
        if ( !approx_equal(K,y1,y2,K*1e-14) )
            N_errors++;
    }
    delete [] A;
    delete [] x;
    delete [] y;
    delete [] y1;
    delete [] y2;
    return N_errors>0;
}

// Test dgemm
static bool test_dgemm( int N )
{
    const int K = TEST_SIZE_MAT;
    double *A = new double[K*K];
    double *B = new double[K*K];
    double *C = new double[K*K];
    double *C1 = new double[K*K];
    double *C2 = new double[K*K];
    random(K*K,A);
    random(K*K,B);
    random(K*K,C);
    const double alpha = 3.141592653589793;
    const double beta  = 1.414213562373095;
    int N_errors = 0;
    for (int i=0; i<K; i++) {
        for (int j=0; j<K; j++) {
            C1[i+j*K] = beta*C[i+j*K];
            for (int k=0; k<K; k++)
                C1[i+j*K] += alpha*A[i+k*K]*B[k+j*K];
        }
    }
    for (int i=0; i<N; i++) {
        memcpy(C2,C,K*K*sizeof(double));
        Lapack::dgemm('N','N',K,K,K,alpha,A,K,B,K,beta,C2,K);
        if ( !approx_equal(K*K,C1,C2,K*1e-14) )
            N_errors++;
    }
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] C1;
    delete [] C2;
    return N_errors>0;
}

// Test dgesv
static bool test_dgesv( int N )
{
    // Test solving a diagonal matrix
    const int K = TEST_SIZE_MAT;
    double *A  = new double[K*K];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    memset(A,0,K*K*sizeof(double));
    random(K,x2);
    random(K,b);
    for (int k=0; k<K; k++) {
        A[k+k*K] = x2[k]+1e-16;
        x1[k] = b[k]/(x2[k]+1e-16);
    }
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(x2,b,K*sizeof(double));
        int error = 0;
        Lapack::dgesv(K,1,A,K,IPIV,x2,K,error);
        N_errors += error==0 ? 0:1;
        if ( !approx_equal(K,x1,x2,K*1e-14) )
            N_errors++;
    }
    delete [] A;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgtsv
static bool test_dgtsv( int N )
{
    // Test solving a tri-diagonal matrix by comparing to dgtsv
    const int K = TEST_SIZE_tri;
    double *A  = new double[K*K];
    double *D  = new double[K];
    double *D2 = new double[K];
    double *DL = new double[K-1];
    double *DL2= new double[K-1];
    double *DU = new double[K-1];
    double *DU2= new double[K-1];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    memset(A,0,K*K*sizeof(double));
    random(K,D);
    random(K-1,DL);
    random(K-1,DU);
    random(K,b);
    A[0] = D[0];
    for (int k=1; k<K; k++) {
        A[k+k*K]     = D[k];
        A[(k-1)+k*K] = DU[k-1];
        A[k+(k-1)*K] = DL[k-1];
    }
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    Lapack::dgesv(K,1,A,K,IPIV,x1,K,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(x2,b,K*sizeof(double));
        memcpy(D2,D,K*sizeof(double));
        memcpy(DL2,DL,(K-1)*sizeof(double));
        memcpy(DU2,DU,(K-1)*sizeof(double));
        Lapack::dgtsv(K,1,DL2,D2,DU2,x2,K,error);
        N_errors += error==0 ? 0:1;
        double err = L2Error(N,x1,x2);
        double norm = L2Norm(N,x1);
        if ( err > 1e-14*norm )
            N_errors++;
    }
    delete [] A;
    delete [] D;
    delete [] D2;
    delete [] DL;
    delete [] DL2;
    delete [] DU;
    delete [] DU2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}
// Test dgtsv
static bool test_dgbsv( int N )
{
    // Test solving a banded-diagonal matrix by comparing to dgtsv
    //    N = 6, KL = 2, KU = 1:
    //        *    *    *    +    +    +
    //        *    *    +    +    +    +
    //        *   a12  a23  a34  a45  a56
    //       a11  a22  a33  a44  a55  a66
    //       a21  a32  a43  a54  a65   *
    //       a31  a42  a53  a64   *    *
    const int K = TEST_SIZE_tri;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2*KL+KU+1;
    double *A  = new double[K*K];
    double *AB = new double[K*K2];
    double *AB2= new double[K*K2];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    random(K*K2,AB);
    random(K,b);
    memset(A,0,K*K*sizeof(double));
    for (int k=0; k<K; k++) {
        for (int k2=-KL; k2<=KU; k2++) {
            if ( k+k2<0 || k+k2>=K ) { continue; }
            A[k+k2+k*K] = AB[k2+2*KL+k*K2];
        }
    }
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    Lapack::dgesv(K,1,A,K,IPIV,x1,K,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(x2,b,K*sizeof(double));
        memcpy(AB2,AB,K*K2*sizeof(double));
        Lapack::dgbsv(K,KL,KU,1,AB2,K2,IPIV,x2,K,error);
        N_errors += error==0 ? 0:1;
        if ( !approx_equal(K,x1,x2,4*K*1e-14) )
            N_errors++;
    }
    delete [] A;
    delete [] AB;
    delete [] AB2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgetrf
static bool test_dgetrf( int N )
{
    // Check dgetrf by performing a factorization and solve and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    double *A  = new double[K*K];
    double *A2 = new double[K*K];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV = new int[K];
    random(K*K,A);
    random(K,b);
    int error = 0;
    memcpy(A2,A,K*K*sizeof(double));
    memcpy(x1,b,K*sizeof(double));
    Lapack::dgesv(K,1,A2,K,IPIV,x1,K,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(A2,A,K*K*sizeof(double));
        Lapack::dgetrf(K,K,A2,K,IPIV,error);
        N_errors += error==0 ? 0:1;
    }
    memcpy(x2,b,K*sizeof(double));
    Lapack::dgetrs('N',K,1,A2,K,IPIV,x2,K,error);
    if ( !approx_equal(K,x1,x2,K*1e-14) )
        N_errors++;
    delete [] A;
    delete [] A2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgttrf
static bool test_dgttrf( int N )
{
    // Check dgttrf by performing a factorization and solve and comparing to dgtsv
    const int K = TEST_SIZE_tri;
    double *D  = new double[K];
    double *D2 = new double[K];
    double *DL = new double[K-1];
    double *DL2= new double[K-1];
    double *DU = new double[K-1];
    double *DU2= new double[K-1];
    double *DU3= new double[K-2];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    random(K,D);
    random(K-1,DL);
    random(K-1,DU);
    random(K,b);
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    memcpy(D2,D,K*sizeof(double));
    memcpy(DL2,DL,(K-1)*sizeof(double));
    memcpy(DU2,DU,(K-1)*sizeof(double));
    Lapack::dgtsv(K,1,DL2,D2,DU2,x1,K,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(D2,D,K*sizeof(double));
        memcpy(DL2,DL,(K-1)*sizeof(double));
        memcpy(DU2,DU,(K-1)*sizeof(double));
        Lapack::dgttrf(K,DL2,D2,DU2,DU3,IPIV,error);
        N_errors += error==0 ? 0:1;
    }
    memcpy(x2,b,K*sizeof(double));
    Lapack::dgttrs('N',K,1,DL2,D2,DU2,DU3,IPIV,x2,K,error);
    if ( !approx_equal(K,x1,x2,K*1e-14) )
        N_errors++;
    delete [] D;
    delete [] D2;
    delete [] DL;
    delete [] DL2;
    delete [] DU;
    delete [] DU2;
    delete [] DU3;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgbtrf
static bool test_dgbtrf( int N )
{
    // Check dgbtrf by performing a factorization and solve and comparing to dgbsv
    const int K = TEST_SIZE_tri;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2*KL+KU+1;
    double *AB = new double[K*K2];
    double *AB2= new double[K*K2];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    random(K*K2,AB);
    random(K,b);
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    memcpy(AB2,AB,K*K2*sizeof(double));
    Lapack::dgbsv(K,KL,KU,1,AB2,K2,IPIV,x1,K,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(AB2,AB,K*K2*sizeof(double));
        Lapack::dgbtrf(K,K,KL,KU,AB2,K2,IPIV,error);
        N_errors += error==0 ? 0:1;
    }
    memcpy(x2,b,K*sizeof(double));
    Lapack::dgbtrs('N',K,KL,KU,1,AB2,K2,IPIV,x2,K,error);
    if ( !approx_equal(K,x1,x2,K*1e-14) )
        N_errors++;
    delete [] AB;
    delete [] AB2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgetrs
static bool test_dgetrs( int N )
{
    // Check dgetrs by performing a factorization and solve and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    double *A  = new double[K*K];
    double *A2 = new double[K*K];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV = new int[K];
    random(K*K,A);
    random(K,b);
    int error = 0;
    memcpy(A2,A,K*K*sizeof(double));
    memcpy(x1,b,K*sizeof(double));
    Lapack::dgesv(K,1,A2,K,IPIV,x1,K,error);
    int N_errors = 0;
    Lapack::dgetrf(K,K,A,K,IPIV,error);
    for (int i=0; i<N; i++) {
        memcpy(A2,A,K*K*sizeof(double));
        memcpy(x2,b,K*sizeof(double));
        Lapack::dgetrs('N',K,1,A2,K,IPIV,x2,K,error);
        N_errors += error==0 ? 0:1;
        if ( !approx_equal(K,x1,x2,K*1e-14) )
            N_errors++;
    }
    delete [] A;
    delete [] A2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgttrs
static bool test_dgttrs( int N )
{
    // Check dgttrs by performing a factorization and solve and comparing to dgtsv
    const int K = TEST_SIZE_tri;
    double *D  = new double[K];
    double *D2 = new double[K];
    double *DL = new double[K-1];
    double *DL2= new double[K-1];
    double *DU = new double[K-1];
    double *DU2= new double[K-1];
    double *DU3= new double[K-2];
    double *DU4= new double[K-2];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    random(K,D);
    random(K-1,DL);
    random(K-1,DU);
    random(K,b);
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    memcpy(D2,D,K*sizeof(double));
    memcpy(DL2,DL,(K-1)*sizeof(double));
    memcpy(DU2,DU,(K-1)*sizeof(double));
    Lapack::dgtsv(K,1,DL2,D2,DU2,x1,K,error);
    Lapack::dgttrf(K,DL,D,DU,DU3,IPIV,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(D2,D,K*sizeof(double));
        memcpy(DL2,DL,(K-1)*sizeof(double));
        memcpy(DU2,DU,(K-1)*sizeof(double));
        memcpy(DU4,DU3,(K-2)*sizeof(double));
        memcpy(x2,b,K*sizeof(double));
        Lapack::dgttrs('N',K,1,DL2,D2,DU2,DU4,IPIV,x2,K,error);
        N_errors += error==0 ? 0:1;
        if ( !approx_equal(K,x1,x2,K*1e-14) )
            N_errors++;
    }
    delete [] D;
    delete [] D2;
    delete [] DL;
    delete [] DL2;
    delete [] DU;
    delete [] DU2;
    delete [] DU3;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgbtrs
static bool test_dgbtrs( int N )
{
    // Check dgbtrs by performing a factorization and solve and comparing to dgbsv
    const int K = TEST_SIZE_tri;
    const int KL = 2;
    const int KU = 2;
    const int K2 = 2*KL+KU+1;
    double *AB = new double[K*K2];
    double *AB2= new double[K*K2];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV  = new int[K];
    random(K*K2,AB);
    random(K,b);
    int error = 0;
    memcpy(x1,b,K*sizeof(double));
    memcpy(AB2,AB,K*K2*sizeof(double));
    Lapack::dgbsv(K,KL,KU,1,AB2,K2,IPIV,x1,K,error);
    Lapack::dgbtrf(K,K,KL,KU,AB,K2,IPIV,error);
    int N_errors = 0;
    for (int i=0; i<N; i++) {
        memcpy(AB2,AB,K*K2*sizeof(double));
        memcpy(x2,b,K*sizeof(double));
        Lapack::dgbtrs('N',K,KL,KU,1,AB2,K2,IPIV,x2,K,error);
        N_errors += error==0 ? 0:1;
        if ( !approx_equal(K,x1,x2,K*1e-14) )
            N_errors++;
    }
    delete [] AB;
    delete [] AB2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    return N_errors>0;
}

// Test dgetri
static bool test_dgetri( int N )
{
    // Check dgetri by performing a factorization, calculating the inverse,
    //   multiplying the rhs, and comparing to dgesv
    const int K = TEST_SIZE_MAT;
    const int LWORK = 8*K;
    double *A  = new double[K*K];
    double *A2 = new double[K*K];
    double *x1 = new double[K];
    double *x2 = new double[K];
    double *b  = new double[K];
    int *IPIV = new int[K];
    double *WORK = new double[LWORK];
    random(K*K,A);
    random(K,b);
    int error = 0;
    memcpy(A2,A,K*K*sizeof(double));
    memcpy(x1,b,K*sizeof(double));
    Lapack::dgesv(K,1,A2,K,IPIV,x1,K,error);
    int N_errors = 0;
    Lapack::dgetrf(K,K,A,K,IPIV,error);
    for (int i=0; i<N; i++) {
        // Compute the inverse
        memcpy(A2,A,K*K*sizeof(double));
        Lapack::dgetri(K,A2,K,IPIV,WORK,LWORK,error);
        N_errors += error==0 ? 0:1;
        // Perform the mat-vec
        memset(x2,0xB6,K*sizeof(double));
        Lapack::dgemv('N',K,K,1,A2,K,b,1,0,x2,1);
        // Check the result
        if ( !approx_equal(K,x1,x2,K*1e-12) )
            N_errors++;
    }
    delete [] A;
    delete [] A2;
    delete [] x1;
    delete [] x2;
    delete [] b;
    delete [] IPIV;
    delete [] WORK;
    return N_errors>0;
}


/******************************************************************
* Some inline functions to acquire/release a mutex                *
******************************************************************/
#ifdef WINDOWS
    HANDLE LapackWrappers_lock_queue = CreateMutex(NULL, FALSE, NULL);
#elif defined(LINUX)
    pthread_mutex_t LapackWrappers_lock_queue;
    int LapackWrappers_lock_queue_error = pthread_mutex_init(&LapackWrappers_lock_queue,NULL);
#else
    #error Not programmed
#endif
#ifdef WINDOWS
    void AMP::Lapack::get_lock( ) {
        WaitForSingleObject( LapackWrappers_lock_queue, INFINITE );
    }
#elif defined(LINUX)
    void AMP::Lapack::get_lock( ) {
        int error = pthread_mutex_lock(&LapackWrappers_lock_queue);
        if ( error == -1 )
            printf("Error locking mutex");
    }
#else
    #error Not programmed
#endif
#ifdef WINDOWS
    void AMP::Lapack::release_lock( ) {
        bool success = ReleaseMutex(LapackWrappers_lock_queue)!=0;
        if ( !success )
            printf("Error unlocking mutex");
    }
#elif defined(LINUX)
    void AMP::Lapack::release_lock( ) {
        int error = pthread_mutex_unlock(&LapackWrappers_lock_queue);
        if ( error == -1 )
            printf("Error unlocking mutex");
    }
#else
    #error Not programmed
#endif



