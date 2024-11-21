#include "AMP/utils/AMPManager.h"
#include "AMP/utils/CopyCast.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"

#include <cstdlib>
#include <memory>
#include <stdio.h>
#ifdef USE_DEVICE
    #include <thrust/device_ptr.h>
#endif

template<typename TypeIn,
         typename TypeOut,
         class AllocatorIn,
         class AllocatorOut,
         typename ptrIn,
         typename ptrOut>
void runTest( AMP::UnitTest &ut, const std::string &pass_msg )
{
    AllocatorIn alloc1;
    TypeIn *v1 = alloc1.allocate( 3 );
    AllocatorOut alloc2;
    TypeOut *v2 = alloc2.allocate( 3 );
    ptrIn v1p( v1 );
    ptrOut v2p( v2 );
    std::srand( 12 );
    for ( int i = 0; i < 3; i++ )
        v1p[i] = static_cast<TypeIn>( std::rand() ) / static_cast<TypeIn>( RAND_MAX );

    // Perform copy-cast
    AMP::Utilities::copyCast<TypeIn, TypeOut>( 3, v1, v2 );
    bool pass = true;
    auto tol  = std::numeric_limits<float>::epsilon();
    for ( int i = 0; i < 3; i++ ) {
        const double err = std::abs( ( v1p[i] - v2p[i] ) / v1p[i] );
        if ( err > tol ) {
            ut.failure( AMP::Utilities::stringf(
                "Cast precision loss %e larger than %e, for entry %d", err, tol, i ) );
            pass = false;
        }
    }
    if ( pass )
        ut.passes( pass_msg );
}


template<class AllocatorIn, class AllocatorOut, typename ptrIn>
void testOverflow( AMP::UnitTest &ut )
{
    AllocatorIn alloc1;
    double *v1 = alloc1.allocate( 3 );
    AllocatorOut alloc2;
    float *v2 = alloc2.allocate( 3 );
    ptrIn v1p( v1 );
    std::srand( 12 );
    for ( int i = 0; i < 3; i++ )
        v1p[i] = static_cast<double>( std::rand() ) / static_cast<double>( RAND_MAX );
    v1p[1] = std::numeric_limits<float>::max() * 2;
    // Perform copy-cast
    try {
        AMP::Utilities::copyCast<double, float>( 3, v1, v2 );
        ut.failure( "CopyCast didn't catch an overflow." );
    } catch ( ... ) {
        ut.passes( "Overflow test succeeded." );
    }
}

int main( int argc, char *argv[] )
{
    // Start AMP (initializes Kokkos)
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    // Test Host
    runTest<double,
            float,
            AMP::HostAllocator<double>,
            AMP::HostAllocator<float>,
            double *,
            float *>( ut, "Host double->float passed" );
    runTest<float,
            double,
            AMP::HostAllocator<float>,
            AMP::HostAllocator<double>,
            float *,
            double *>( ut, "Host float->double passed" );

#ifndef USE_DEVICE
    // Test Overflow
    testOverflow<AMP::HostAllocator<double>, AMP::HostAllocator<float>, double *>( ut );
#endif

#ifdef USE_DEVICE
    // Test Managed
    runTest<double,
            float,
            AMP::ManagedAllocator<double>,
            AMP::ManagedAllocator<float>,
            thrust::device_ptr<double>,
            thrust::device_ptr<float>>( ut, "Managed double->float passed" );
    runTest<float,
            double,
            AMP::ManagedAllocator<float>,
            AMP::ManagedAllocator<double>,
            thrust::device_ptr<float>,
            thrust::device_ptr<double>>( ut, "Managed float->double passed" );
    // Test Managed
    runTest<double,
            float,
            AMP::DeviceAllocator<double>,
            AMP::DeviceAllocator<float>,
            thrust::device_ptr<double>,
            thrust::device_ptr<float>>( ut, "Device double->float passed" );
    runTest<float,
            double,
            AMP::DeviceAllocator<float>,
            AMP::DeviceAllocator<double>,
            thrust::device_ptr<float>,
            thrust::device_ptr<double>>( ut, "Device float->double passed" );
#endif

    int N_errors = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_errors;
}
