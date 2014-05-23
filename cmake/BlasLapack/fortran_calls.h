#ifndef FORTRAN_LAPACK_CALLS
#define FORTRAN_LAPACK_CALLS

extern "C" {

#if defined(_WIN32) || defined(__hpux)
#define FORTRAN_WRAPPER(x) x
#else
#define FORTRAN_WRAPPER(x) x ## _
#endif

// Misc
extern void FORTRAN_WRAPPER(dfill)( int*, double*, double* );

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

};

#endif


