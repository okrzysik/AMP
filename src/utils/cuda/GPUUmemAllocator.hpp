#ifndef included_GPUUmemAllocator_HPP_
#define included_GPUUmemAllocator_HPP_


#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

#include "utils/Utilities.h"

inline void checkError(cudaError_t e){
    if(e != cudaSuccess){
        std::cout<<cudaGetErrorString(e) <<std::endl;
        exit(e);
    }
}

template <typename T>
T* GPUUmemAllocator<T>::allocate(size_t n)
{
    T* d_m;
    cudaError_t E = cudaMallocManaged(&d_m,sizeof(T)*n);
    checkError(E);  
    
    return d_m;
}

template <typename T>
void GPUUmemAllocator<T>::deallocate(T* p,size_t n)
{
    (void)n;
    cudaError_t E = cudaFree(p);
    checkError(E);  
} 

template <typename T>
void GPUUmemAllocator<T>::construct(T* p)
{
    new(p) T;
}

template <typename T>
void GPUUmemAllocator<T>::destroy(T* p)
{
    NULL_USE(p);
//    delete(p);
}

#endif
