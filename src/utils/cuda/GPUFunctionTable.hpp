#ifndef included_AMP_GPUFunctionTable_HPP_
#define included_AMP_GPUFunctionTable_HPP_

#include "cublas_v2.h"
#include "curand.h"


#include "AMP/utils/Array.h"
#include "AMP/utils/UtilityMacros.h"


namespace AMP {

// Kernel Wrappers
template<class TYPE>
void transformReLUW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
void transformAbsW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
void transformTanhW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
void transformHardTanhW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
void transformSigmoidW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
void transformSoftPlusW( const TYPE *d_a, TYPE *d_b, size_t n );

template<class TYPE>
TYPE sumW( const TYPE *d_a, size_t n );

template<class TYPE>
bool equalsW( const TYPE *d_a, const TYPE *d_b, TYPE tol, size_t n );

// Rand functions
template<class TYPE, class FUN, class ALLOC>
inline void GPUFunctionTable::rand( Array<TYPE, FUN, ALLOC> &x )
{
    rand<TYPE>( x.length(), x.data() );
}

template<>
inline void GPUFunctionTable::rand<int>( size_t n, int *d_x )
{
    curandGenerator_t gen;
    curandCreateGenerator( &gen, CURAND_RNG_PSEUDO_DEFAULT );
    curandSetPseudoRandomGeneratorSeed( gen, time( NULL ) );
    curandGenerate( gen, (unsigned int *) d_x, n );
    curandDestroyGenerator( gen );
}

template<>
inline void GPUFunctionTable::rand<float>( size_t n, float *d_x )
{
    curandGenerator_t gen;
    curandCreateGenerator( &gen, CURAND_RNG_PSEUDO_DEFAULT );
    curandSetPseudoRandomGeneratorSeed( gen, time( NULL ) );
    curandGenerateUniform( gen, d_x, n );
    curandDestroyGenerator( gen );
}

template<>
inline void GPUFunctionTable::rand<double>( size_t n, double *d_x )
{
    curandGenerator_t gen;
    curandCreateGenerator( &gen, CURAND_RNG_PSEUDO_DEFAULT );
    curandSetPseudoRandomGeneratorSeed( gen, time( NULL ) );
    curandGenerateUniformDouble( gen, d_x, n );
    curandDestroyGenerator( gen );
}


// Specialized transform functions - temporary solution
template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformReLU( const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformReLUW<TYPE>( A.data(), B.data(), A.length() );
}

template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformAbs( const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformAbsW<TYPE>( A.data(), B.data(), A.length() );
}
template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformTanh( const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformTanhW<TYPE>( A.data(), B.data(), A.length() );
}

template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformHardTanh( const Array<TYPE, FUN, ALLOC> &A,
                                          Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformHardTanhW<TYPE>( A.data(), B.data(), A.length() );
}

template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformSigmoid( const Array<TYPE, FUN, ALLOC> &A,
                                         Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformSigmoidW<TYPE>( A.data(), B.data(), A.length() );
}

template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::transformSoftPlus( const Array<TYPE, FUN, ALLOC> &A,
                                          Array<TYPE, FUN, ALLOC> &B )
{
    B.resize( A.size() );
    transformSoftPlusW<TYPE>( A.data(), B.data(), A.length() );
}

// Specialized reductions
template<class TYPE, class FUN, class ALLOC>
TYPE GPUFunctionTable::sum( const Array<TYPE, FUN, ALLOC> &A )
{
    if ( A.length() == 0 ) {
        return TYPE();
    }
    return sumW<TYPE>( A.data(), A.length() );
}

template<class TYPE, class FUN, class ALLOC>
bool GPUFunctionTable::equals( const Array<TYPE, FUN, ALLOC> &A,
                               const Array<TYPE, FUN, ALLOC> &B,
                               TYPE tol )
{
    bool eq = true;
    AMP_INSIST( A.size() == B.size(), "Sizes of A and B do not match" );
    eq = equalsW( A.data(), B.data(), tol, A.length() );

    return eq;
}

// GEMM wrappers - might not need these
template<class TYPE>
inline void gemmWrapper( char TRANSA,
                         char TRANSB,
                         int M,
                         int N,
                         int K,
                         TYPE alpha,
                         const TYPE *A,
                         int LDA,
                         const TYPE *B,
                         int LDB,
                         TYPE beta,
                         TYPE *C,
                         int LDC )
{
    AMP_ERROR( "Not supported for the provided type" );
}

template<>
inline void gemmWrapper<float>( char TRANSA,
                                char TRANSB,
                                int M,
                                int N,
                                int K,
                                float alpha,
                                const float *A,
                                int LDA,
                                const float *B,
                                int LDB,
                                float beta,
                                float *C,
                                int LDC )
{
    cublasOperation_t transa, transb;
    if ( TRANSA == 'N' ) {
        transa = CUBLAS_OP_N;
    } else if ( TRANSA == 'T' ) {
        transa = CUBLAS_OP_T;
    } else {
        AMP_ERROR( "Invalid specification" );
    }
    if ( TRANSB == 'N' ) {
        transb = CUBLAS_OP_N;
    } else if ( TRANSB == 'T' ) {
        transb = CUBLAS_OP_T;
    } else {
        AMP_ERROR( "Invalid specification" );
    }

    cublasHandle_t handle;
    cublasCreate( &handle );
    cublasSgemm( handle, transa, transb, M, N, K, &alpha, A, LDA, B, LDB, &beta, C, LDC );
}

template<>
inline void gemmWrapper<double>( char TRANSA,
                                 char TRANSB,
                                 int M,
                                 int N,
                                 int K,
                                 double alpha,
                                 const double *A,
                                 int LDA,
                                 const double *B,
                                 int LDB,
                                 double beta,
                                 double *C,
                                 int LDC )
{

    cublasOperation_t transa, transb;
    if ( TRANSA == 'N' ) {
        transa = CUBLAS_OP_N;
    } else if ( TRANSA == 'T' ) {
        transa = CUBLAS_OP_T;
    } else {
        AMP_ERROR( "Invalid specification" );
    }
    if ( TRANSB == 'N' ) {
        transb = CUBLAS_OP_N;
    } else if ( TRANSB == 'T' ) {
        transb = CUBLAS_OP_T;
    } else {
        AMP_ERROR( "Invalid specification" );
    }

    cublasHandle_t handle;
    cublasCreate( &handle );
    cublasDgemm( handle, transa, transb, M, N, K, &alpha, A, LDA, B, LDB, &beta, C, LDC );
}


/* Functions not yet implemented */

template<class TYPE, class FUN, class ALLOC, typename LAMBDA>
inline void
GPUFunctionTable::transform( LAMBDA &, const Array<TYPE, FUN, ALLOC> &, Array<TYPE, FUN, ALLOC> & )
{
    AMP_ERROR( "Not implemented for GPU" );
}

template<class TYPE, class FUN, class ALLOC, typename LAMBDA>
inline void GPUFunctionTable::transform( LAMBDA &,
                                         const Array<TYPE, FUN, ALLOC> &,
                                         const Array<TYPE, FUN, ALLOC> &,
                                         Array<TYPE, FUN, ALLOC> & )
{
    AMP_ERROR( "Not implemented for GPU" );
}

template<class TYPE, class FUN, class ALLOC, typename LAMBDA>
inline TYPE GPUFunctionTable::reduce( LAMBDA &, const Array<TYPE, FUN, ALLOC> &, const TYPE & )
{
    AMP_ERROR( "Not implemented for GPU" );
    return 0;
}

template<class TYPE, class FUN, class ALLOC, typename LAMBDA>
inline TYPE GPUFunctionTable::reduce( LAMBDA &,
                                      const Array<TYPE, FUN, ALLOC> &,
                                      const Array<TYPE, FUN, ALLOC> &,
                                      const TYPE & )
{
    AMP_ERROR( "Not implemented for GPU" );
    return 0;
}

template<class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::multiply( const Array<TYPE, FUN, ALLOC> &,
                                 const Array<TYPE, FUN, ALLOC> &,
                                 Array<TYPE, FUN, ALLOC> & )
{
    AMP_ERROR( "not implemented" );
}
} // namespace AMP
#endif
