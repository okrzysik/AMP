#ifndef FORTRAN_LAPACK_CALLS
#define FORTRAN_LAPACK_CALLS

extern "C" {

#if defined(_WIN32) || defined(__hpux)
#define FORTRAN_WRAPPER(x) x
#else
#define FORTRAN_WRAPPER(x) x ## _
#endif

// double precision functions
// Misc
extern void   FORTRAN_WRAPPER( dfill)( int*, double*, double* );
extern int    FORTRAN_WRAPPER(idamax)( int* N, double* X, int* INCX );
extern double FORTRAN_WRAPPER(dlamch)( char* cmach );

// Level 1 BLAS Routines
extern void   FORTRAN_WRAPPER(dswap)( int* N,                double* x, int* INCX, double* y, int* INCY );   //  x <-> y
extern void   FORTRAN_WRAPPER(dscal)( int* N, double* alpha, double* x, int* INCX );                         //  x = alpha*x
extern void   FORTRAN_WRAPPER(dcopy)( int* N,                double* x, int* INCX, double* y, int* INCY );   //  y = x
extern void   FORTRAN_WRAPPER(daxpy)( int* N, double* alpha, double* x, int* INCX, double* y, int* INCY );   //  y = alpha*x + y
extern double FORTRAN_WRAPPER( ddot)( int* N,                double* x, int* INCX, double* y, int* INCY );   //  dot = xT*y
extern double FORTRAN_WRAPPER(dasum)( int* N,                double* x, int* INCX                       );   //  asum = sum(|x|)
extern double FORTRAN_WRAPPER(damax)( int* N,                double* x, int* INCX                       );   //  amax = max(|x|)
extern double FORTRAN_WRAPPER(dnrm2)( int* N, double* X, int* INCX );                                        // dnrm2 = sqrt(x*x)

// Level 2 BLAS Routines 
extern void FORTRAN_WRAPPER(dgemv)( char* TRANS, int* M, int* N, double* alpha, double* A, int* LDA, double* x, int* INCX, double* beta, double* y, int* INCY );  //  y = alpha*A*x + beta*y
extern void FORTRAN_WRAPPER(zgemv)( char* TRANS, int* M, int* N, double* alpha, double* A, int* LSA, double* x, int* INCX, double* beta, double* y, int* INCY );  //  y = alpha*A*x + beta*y
extern void FORTRAN_WRAPPER( dger)( int* M, int* N, double* alpha, double* x, int* INCX, double* y, int* INCY, double* A, int* LDA );                             //  A = alpha*x*yT + A
extern void FORTRAN_WRAPPER(zgeru)( int* M, int* N, double* alpha, double* x, int* INCX, double* y, int* INCY, double* A, int* LDA );                             //  A = alpha*x*yT + A

// Level 3 BLAS Routines 
extern void FORTRAN_WRAPPER(dgemm)( char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha, double*A, int* LDA, double* B, int* LDB, double* beta, double* C, int* LDC );  //  C = alpha*op(A)*op(B) + beta*C
extern void FORTRAN_WRAPPER(zgemm)( char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* alpha, double*A, int* LDA, double* B, int* LDB, double* beta, double* C, int* LDC );  //  C = alpha*op(A)*op(B) + beta*C

// LAPACK Routines 
extern void FORTRAN_WRAPPER( dgesv)( int* N, int* NRHS, double* A, int* LDA , int* IPIV, double* B, int* LDB, int* info );                                        //  Solve Ax=b using LU decompositionon for a general M-by-N system
extern void FORTRAN_WRAPPER( dgtsv)( int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO );                                      //  Solve Ax=b using LU decompositionon for a tridiagonal system
extern void FORTRAN_WRAPPER( dgbsv)( int* N, int *KL, int *KU, int* NRHS, double* AB, int* LDA , int* IPIV, double* B, int* LDB, int* info );                     //  Solve Ax=b using LU decompositionon for a band matrix
extern void FORTRAN_WRAPPER( dgbmv)( char* TRANS, int* M, int* N, int* KL, int* KU, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
extern void FORTRAN_WRAPPER(dgetrf)( int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO );                                                                 //  Compute LU factorization on a general M-by-N system
extern void FORTRAN_WRAPPER(dgttrf)( int* N, double* DL, double* D, double* DU, double* DU2, int* IPIV, int* INFO );                                              //  Compute LU factorization on a tridiagonal system
extern void FORTRAN_WRAPPER(dgbtrf)( int* M, int *N, int *KL, int *KU, double *AB, int *LDAB, int* IPIV, int* INFO );                                             //  Compute LU factorization on a band matrix  
extern void FORTRAN_WRAPPER(dgetrs)( char* TRANS, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO );                            //  Solve Ax=b using the LU factorization computed by dgetrf for a general M-by-N system
extern void FORTRAN_WRAPPER(dgttrs)( char* TRANS, int* N, int* NRHS, double* DL, double* D, double* DU, double* DU2, int* IPIV, double* B, int* LDB, int* INFO ); //  Solve Ax=b using the LU factorization computed by dgttrf for a tridiagonal system
extern void FORTRAN_WRAPPER(dgbtrs)( char* TRANS, int* N, int *KL, int *KU, int* NRHS, double* AB, int* LDAB, int* IPIV, double* B, int* LDB, int* INFO );        //  Solve Ax=b using the LU factorization computed by dgbtrf for a band matrix
extern void FORTRAN_WRAPPER(dgetri)( int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO );                                               //  Compute the inverse matrix using the LU factorization computed by dgetrf 
extern void FORTRAN_WRAPPER(dtrsm)( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, double* alpha, double* A, int* LDA, double* B, int* LDB );  //  Solve one of the matrix equations op(A)*X=alpha*B or X*op(A)=alpha*B

// single precision functions
// Misc
extern void   FORTRAN_WRAPPER( sfill)( int*, float*, float* );
extern int    FORTRAN_WRAPPER(isamax)( int* N, float* X, int* INCX );
extern float FORTRAN_WRAPPER( slamch)( char* cmach );

// Level 1 BLAS Routines
extern void   FORTRAN_WRAPPER(sswap)( int* N,                float* x, int* INCX, float* y, int* INCY );   //  x <-> y
extern void   FORTRAN_WRAPPER(sscal)( int* N, float* alpha, float* x, int* INCX );                         //  x = alpha*x
extern void   FORTRAN_WRAPPER(scopy)( int* N,                float* x, int* INCX, float* y, int* INCY );   //  y = x
extern void   FORTRAN_WRAPPER(saxpy)( int* N, float* alpha, float* x, int* INCX, float* y, int* INCY );   //  y = alpha*x + y
extern float FORTRAN_WRAPPER( sdot)( int* N,                float* x, int* INCX, float* y, int* INCY );   //  dot = xT*y
extern float FORTRAN_WRAPPER(sasum)( int* N,                float* x, int* INCX                       );   //  asum = sum(|x|)
extern float FORTRAN_WRAPPER(samax)( int* N,                float* x, int* INCX                       );   //  amax = max(|x|)
extern float FORTRAN_WRAPPER(snrm2)( int* N, float* X, int* INCX );                                        // dnrm2 = sqrt(x*x)

// Level 2 BLAS Routines 
extern void FORTRAN_WRAPPER(sgemv)( char* TRANS, int* M, int* N, float* alpha, float* A, int* LDA, float* x, int* INCX, float* beta, float* y, int* INCY );  //  y = alpha*A*x + beta*y
    //extern void FORTRAN_WRAPPER(zgemv)( char* TRANS, int* M, int* N, float* alpha, float* A, int* LSA, float* x, int* INCX, float* beta, float* y, int* INCY );  //  y = alpha*A*x + beta*y
extern void FORTRAN_WRAPPER( sger)( int* M, int* N, float* alpha, float* x, int* INCX, float* y, int* INCY, float* A, int* LDA );                             //  A = alpha*x*yT + A
    //extern void FORTRAN_WRAPPER(zgeru)( int* M, int* N, float* alpha, float* x, int* INCX, float* y, int* INCY, float* A, int* LDA );                             //  A = alpha*x*yT + A

// Level 3 BLAS Routines 
extern void FORTRAN_WRAPPER(sgemm)( char* TRANSA, char* TRANSB, int* M, int* N, int* K, float* alpha, float*A, int* LDA, float* B, int* LDB, float* beta, float* C, int* LDC );  //  C = alpha*op(A)*op(B) + beta*C
    //extern void FORTRAN_WRAPPER(zgemm)( char* TRANSA, char* TRANSB, int* M, int* N, int* K, float* alpha, float*A, int* LDA, float* B, int* LDB, float* beta, float* C, int* LDC );  //  C = alpha*op(A)*op(B) + beta*C

// LAPACK Routines 
extern void FORTRAN_WRAPPER( sgesv)( int* N, int* NRHS, float* A, int* LDA , int* IPIV, float* B, int* LDB, int* info );                                        //  Solve Ax=b using LU decompositionon for a general M-by-N system
extern void FORTRAN_WRAPPER( sgtsv)( int* N, int* NRHS, float* DL, float* D, float* DU, float* B, int* LDB, int* INFO );                                      //  Solve Ax=b using LU decompositionon for a tridiagonal system
extern void FORTRAN_WRAPPER( sgbsv)( int* N, int *KL, int *KU, int* NRHS, float* AB, int* LDA , int* IPIV, float* B, int* LDB, int* info );                     //  Solve Ax=b using LU decompositionon for a band matrix
extern void FORTRAN_WRAPPER( sgbmv)( char* TRANS, int* M, int* N, int* KL, int* KU, float* ALPHA, float* A, int* LDA, float* X, int* INCX, float* BETA, float* Y, int* INCY);
extern void FORTRAN_WRAPPER(sgetrf)( int* M, int* N, float* A, int* LDA, int* IPIV, int* INFO );                                                                 //  Compute LU factorization on a general M-by-N system
extern void FORTRAN_WRAPPER(sgttrf)( int* N, float* DL, float* D, float* DU, float* DU2, int* IPIV, int* INFO );                                              //  Compute LU factorization on a tridiagonal system
extern void FORTRAN_WRAPPER(sgbtrf)( int* M, int *N, int *KL, int *KU, float *AB, int *LDAB, int* IPIV, int* INFO );                                             //  Compute LU factorization on a band matrix  
extern void FORTRAN_WRAPPER(sgetrs)( char* TRANS, int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO );                            //  Solve Ax=b using the LU factorization computed by dgetrf for a general M-by-N system
extern void FORTRAN_WRAPPER(sgttrs)( char* TRANS, int* N, int* NRHS, float* DL, float* D, float* DU, float* DU2, int* IPIV, float* B, int* LDB, int* INFO ); //  Solve Ax=b using the LU factorization computed by dgttrf for a tridiagonal system
extern void FORTRAN_WRAPPER(sgbtrs)( char* TRANS, int* N, int *KL, int *KU, int* NRHS, float* AB, int* LDAB, int* IPIV, float* B, int* LDB, int* INFO );        //  Solve Ax=b using the LU factorization computed by dgbtrf for a band matrix
extern void FORTRAN_WRAPPER(sgetri)( int* N, float* A, int* LDA, int* IPIV, float* WORK, int* LWORK, int* INFO );                                               //  Compute the inverse matrix using the LU factorization computed by dgetrf 
extern void FORTRAN_WRAPPER(strsm)( char* SIDE, char* UPLO, char* TRANSA, char* DIAG, int* M, int* N, float* alpha, float* A, int* LDA, float* B, int* LDB );  //  Solve one of the matrix equations op(A)*X=alpha*B or X*op(A)=alpha*B

};

#endif


