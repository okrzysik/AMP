#include "AMP/utils/cuda/helper_cuda.h"
#include <iostream>


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

    return 0;
}
