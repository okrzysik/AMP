#include "AMP/utils/cuda/CudaAllocator.h"
#include "AMP/utils/cuda/helper_cuda.h"

#include <iostream>
#include <memory>


int main()
{
    // cudaError_t
    std::cout << "cudaSuccess: " << cudaGetName( cudaSuccess ) << std::endl;
    std::cout << "cudaErrorNoDevice: " << cudaGetName( cudaErrorNoDevice ) << std::endl;
    std::cout << "cudaErrorMemoryAllocation: " << cudaGetName( cudaErrorMemoryAllocation )
              << std::endl;
    std::cout << "cudaErrorMissingConfiguration: " << cudaGetName( cudaErrorMissingConfiguration )
              << std::endl;

    // CUresult
    std::cout << "CUDA_SUCCESS: " << cudaGetName( CUDA_SUCCESS ) << std::endl;
    std::cout << "CUDA_ERROR_NO_DEVICE: " << cudaGetName( CUDA_ERROR_NO_DEVICE ) << std::endl;
    std::cout << "CUDA_ERROR_OUT_OF_MEMORY: " << cudaGetName( CUDA_ERROR_OUT_OF_MEMORY )
              << std::endl;

    // Memory pointer type
    AMP::CudaDevAllocator<double> devAllocator;
    AMP::CudaManagedAllocator<double> managedAllocator;
    AMP::CudaHostAllocator<double> hostAllocator;
    std::allocator<double> stdAllocator;
    size_t N     = 100;
    auto device  = devAllocator.allocate( N );
    auto managed = managedAllocator.allocate( N );
    auto host    = hostAllocator.allocate( N );
    auto std     = stdAllocator.allocate( N );
    std::cout << std::endl;
    std::cout << "Device: " << getString( getMemoryType( device ) ) << std::endl;
    std::cout << "Managed: " << getString( getMemoryType( managed ) ) << std::endl;
    std::cout << "Host: " << getString( getMemoryType( host ) ) << std::endl;
    std::cout << "std: " << getString( getMemoryType( std ) ) << std::endl;
    devAllocator.deallocate( device, N );
    managedAllocator.deallocate( managed, N );
    hostAllocator.deallocate( host, N );
    stdAllocator.deallocate( std, N );
    return 0;
}
