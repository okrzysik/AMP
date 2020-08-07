#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/PIO.h"
#include "AMP/vectors/SimpleVector.h"

#ifdef USE_OPENMP
#include "AMP/vectors/operations/OpenMP/VectorOperationsOpenMP.h"
#endif
#ifdef USE_CUDA
#include "AMP/vectors/data/cuda/VectorDataGPU.h"
#include "AMP/vectors/operations/cuda/VectorOperationsCuda.h"
#endif

#include "ProfilerApp.h"

#include <chrono>


static inline double speedup( size_t x, size_t y )
{
    return static_cast<double>( y ) / static_cast<double>( x );
}

struct test_times {
    size_t clone;
    size_t zero;
    size_t setToScalar;
    size_t setRandomValues;
    size_t L1Norm;
    size_t L2Norm;
    size_t maxNorm;
    size_t axpy;
    size_t dot;
    size_t min;
    size_t max;
    size_t multiply;
    size_t divide;
    void print()
    {
        AMP::pout << "    clone: " << clone / 1000 << " us" << std::endl;
        AMP::pout << "    zero: " << zero / 1000 << " us" << std::endl;
        AMP::pout << "    setToScalar: " << setToScalar / 1000 << " us" << std::endl;
        AMP::pout << "    setRandom: " << setRandomValues / 1000 << " us" << std::endl;
        AMP::pout << "    L1Norm: " << L1Norm / 1000 << " us" << std::endl;
        AMP::pout << "    L2Norm: " << L2Norm / 1000 << " us" << std::endl;
        AMP::pout << "    maxNorm: " << maxNorm / 1000 << " us" << std::endl;
        AMP::pout << "    axpy: " << axpy / 1000 << " us" << std::endl;
        AMP::pout << "    min: " << min / 1000 << " us" << std::endl;
        AMP::pout << "    max: " << max / 1000 << " us" << std::endl;
        AMP::pout << "    dot: " << dot / 1000 << " us" << std::endl;
        AMP::pout << "    multiply: " << multiply / 1000 << " us" << std::endl;
        AMP::pout << "    divide: " << divide / 1000 << " us" << std::endl;
    }
    void print_speedup( const test_times &time0 )
    {
        AMP::pout << "  Speedup: " << std::endl;
        AMP::pout << "    clone: " << speedup( clone, time0.clone ) << std::endl;
        AMP::pout << "    zero: " << speedup( zero, time0.zero ) << std::endl;
        AMP::pout << "    setToScalar: " << speedup( setToScalar, time0.setToScalar ) << std::endl;
        AMP::pout << "    setRandom: " << speedup( setRandomValues, time0.setRandomValues )
                  << std::endl;
        AMP::pout << "    L1Norm: " << speedup( L1Norm, time0.L1Norm ) << std::endl;
        AMP::pout << "    L2Norm: " << speedup( L2Norm, time0.L2Norm ) << std::endl;
        AMP::pout << "    maxNorm: " << speedup( maxNorm, time0.maxNorm ) << std::endl;
        AMP::pout << "    axpy: " << speedup( axpy, time0.axpy ) << std::endl;
        AMP::pout << "    min: " << speedup( min, time0.min ) << std::endl;
        AMP::pout << "    max: " << speedup( max, time0.max ) << std::endl;
        AMP::pout << "    dot: " << speedup( dot, time0.dot ) << std::endl;
        AMP::pout << "    multiply: " << speedup( multiply, time0.multiply ) << std::endl;
        AMP::pout << "    divide: " << speedup( divide, time0.divide ) << std::endl;
    }
};


#define runTest0( TEST )                                                                      \
    do {                                                                                      \
        auto t1 = std::chrono::steady_clock::now();                                           \
        vec->TEST();                                                                          \
        auto t2    = std::chrono::steady_clock::now();                                        \
        times.TEST = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count(); \
    } while ( 0 )
#define runTest1( TEST, X )						\
    do {                                                                                      \
        auto t1 = std::chrono::steady_clock::now();                                           \
        vec->TEST( X );						\
        auto t2    = std::chrono::steady_clock::now();                                        \
        times.TEST = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count(); \
    } while ( 0 )
#define runTest2( TEST, X, Y )                                                                \
    do {                                                                                      \
        auto t1 = std::chrono::steady_clock::now();                                           \
        vec->TEST( X, Y );                                                                    \
        auto t2    = std::chrono::steady_clock::now();                                        \
        times.TEST = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count(); \
    } while ( 0 )
#define runTest3( TEST, X, Y, Z )                                                             \
    do {                                                                                      \
        auto t1 = std::chrono::steady_clock::now();                                           \
        vec->TEST( X, Y, Z );                                                                 \
        auto t2    = std::chrono::steady_clock::now();                                        \
        times.TEST = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count(); \
    } while ( 0 )
#define runTest4( TEST, X, Y, Z, V )					\
    do {                                                                                      \
        auto t1 = std::chrono::steady_clock::now();                                           \
        vec->TEST( X, Y, Z, V );						\
        auto t2    = std::chrono::steady_clock::now();                                        \
        times.TEST = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count(); \
    } while ( 0 )


test_times testPerformance( AMP::LinearAlgebra::Vector::shared_ptr vec )
{
    test_times times;
    // Test the performance of clone
    auto t1     = std::chrono::steady_clock::now();
    auto vec2   = vec->cloneVector();
    auto t2     = std::chrono::steady_clock::now();
    times.clone = std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count();
    auto vec3   = vec->cloneVector();
    vec2->setRandomValues(vec2);
    vec3->setRandomValues(vec3);
    // Run the tests
    runTest1( zero, vec );
    runTest2( setToScalar, 3.14, vec );
    runTest1( setRandomValues, vec );
    runTest1( L1Norm, vec );
    runTest1( L2Norm, vec );
    runTest1( maxNorm, vec );
    runTest1( min, vec );
    runTest1( max, vec );
    runTest2( dot, *vec2, *vec );
    runTest4( axpy, 2.5, *vec2, *vec3, *vec );
    runTest3( multiply, *vec2, *vec3, *vec );
    runTest3( multiply, *vec2, *vec3, *vec );

    return times;
}


int main( int argc, char **argv )
{

    AMP::AMPManager::startup( argc, argv );
    PROFILE_DISABLE();

    {
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        int rank = globalComm.getRank();

        size_t N = 4e6;
        auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "vec" );

        auto vec   = AMP::LinearAlgebra::SimpleVector<double>::create( N, var, globalComm );
        auto time0 = testPerformance( vec );
        if ( rank == 0 ) {
            AMP::pout << "SimpleVector:" << std::endl;
            time0.print();
            AMP::pout << std::endl;
        }

#ifdef USE_OPENMP
        vec = AMP::LinearAlgebra::
            SimpleVector<double, AMP::LinearAlgebra::VectorOperationsOpenMP<double>>::create(
                N, var, globalComm );
        auto time_openmp = testPerformance( vec );
        if ( rank == 0 ) {
            AMP::pout << "SimpleVector<OpenMP>:" << std::endl;
            time_openmp.print();
            time_openmp.print_speedup( time0 );
            AMP::pout << std::endl;
        }
#endif

#ifdef USE_CUDA
        vec = AMP::LinearAlgebra::SimpleVector<
            double,
            AMP::LinearAlgebra::VectorOperationsCuda<double>,
            AMP::LinearAlgebra::VectorDataGPU<double>>::create( N, var, globalComm );
        auto time_cuda = testPerformance( vec );
        if ( rank == 0 ) {
            AMP::pout << "SimpleVector<CUDA>:" << std::endl;
            time_cuda.print();
            time_cuda.print_speedup( time0 );
            AMP::pout << std::endl;
        }
#endif
    }

    AMP::AMPManager::shutdown();
    return 0;
}
