#ifndef included_GPUDevAllocator_HPP_
#define included_GPUDevAllocator_HPP_

#include "utils/cuda/GPUDevAllocator.h"
#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

inline void checkErrorD(cudaError_t e){
    if(e != cudaSuccess){
        std::cerr<<cudaGetErrorString(e) <<std::endl;
        exit(-1);   
    }
}


template <typename T>
T* GPUDevAllocator<T>::allocate(size_t n)
{
    T* d_m;
    cudaError_t E = cudaMalloc(&d_m,sizeof(T)*n);
    checkErrorD(E);  
    
    return d_m;
}

template <typename T>
void GPUDevAllocator<T>::deallocate(T* p,size_t n)
{
    (void)n;
    cudaError_t E = cudaFree(p);
    checkErrorD(E);  

} 

template <typename T>
void GPUDevAllocator<T>::construct(T* p)
{
    new(p) T;
}

template <typename T>
void GPUDevAllocator<T>::destroy(T* p)
{
//    delete(p);
}

#endif
