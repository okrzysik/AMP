#include "AMP/utils/cuda/testGPUAllocators.hpp"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Array.h"
#include "AMP/utils/Array.hpp"
#include "AMP/utils/FunctionTable.h"
#include "AMP/utils/FunctionTable.hpp"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/cuda/CudaAllocator.h"
#include <cuda.h>


int main( int argc, char *argv[] )
{
    // Declare Arrays and utils
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    AMP::Array<double, AMP::FunctionTable, AMP::CudaManagedAllocator<double>> A;
    AMP::Array<double, AMP::FunctionTable, AMP::CudaDevAllocator<double>> B;
    AMP::Array<double> R;

    const size_t n             = 10;
    std::vector<size_t> v1     = { n };
    std::vector<size_t> v2     = { n, n };
    std::vector<size_t> v3     = { n, n, n };
    std::vector<size_t> v4     = { n, n, n, n };
    std::vector<size_t> vEmpty = { 0, 0, 0, 0 };
    A.resize( v1 );
    A.resize( v4 );
    A.resize( v3 );
    A.resize( vEmpty );
    A.resize( v2 );

    // Test to make sure the data is operable
    KernelWrapper<double> K;
    K.setData( A.data(), 4.0, A.length() );
    K.opData( A.data(), A.length() );
    cudaDeviceSynchronize();

    B.resize( v2 );
    K.setData( B.data(), 4.0, B.length() );
    K.opData( B.data(), B.length() );
    R.resize( v2 );
    cudaMemcpy( R.data(), B.data(), sizeof( double ) * R.length(), cudaMemcpyDeviceToHost );
    bool pass = true;
    for ( size_t i = 0; i < R.length(); i++ ) {
        if ( A.data()[i] != R.data()[i] ) {
            pass = false;
        }
    }

    if ( pass ) {
        ut.passes( "GPU Allocator Success" );
    } else {
        ut.failure( "GPU Allocator Failure" );
    }

    ut.report();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
