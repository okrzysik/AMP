#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/FunctionTable.h"
#include "utils/cuda/GPUUmemAllocator.h"
#include "utils/cuda/GPUDevAllocator.h"
#include "utils/cuda/testGPUAllocators.hpp"
#include "utils/Array.h"

/*
template <typename T, class FTable, template <T> class Alloc>
void testGpuAllocators(AMP::UnitTest *ut, AMP::Array<T,FTable,Alloc<T> >& A)
{
    const size_t n = 10;
    const size_t l = (size_t) - 1;
    std::vector<size_t> v1 = {n};
    std::vector<size_t> v2 = {n,n};
    std::vector<size_t> v3 = {n,n,n};
    std::vector<size_t> v4 = {n,n,n,n};
    std::vector<size_t> vLarge = {l,l,l,l};
    std::vector<size_t> vEmpty = {0,0,0,0};
    //allocate a simple array
    A.allocate(v1);
    A.resize(v4); 
    A.resize(v3);
    A.resize(v2);
    //simple type
    //custom type
    //1d 2d 3d and 4d
    //deallocate afterwords
    //do something to it
    //allocate an array thats too large
    A.resize(vLarge);
    A.resize(vEmpty);
    //allocate 0 size array
    //clear the memory
}*/

int main(int argc, char* argv[])
{
    AMP::AMPManager::startup(argc,argv);
    AMP::UnitTest ut; 
    AMP::Array<double, AMP::FunctionTable, GPUUmemAllocator<double> > A;
    AMP::Array<double, AMP::FunctionTable, GPUDevAllocator<double> > B;
  //  AMP::Array<float, AMP::FunctionTable, Gpu_Umem_Allocator<float> > C;
    //AMP::Array<float, AMP::FunctionTable, Gpu_dev_Allocator<float> > D;
    //testGpuAllocators(&ut,A); 
    //testGpuAllocators(&ut,B); 
    //testGpuAllocators(&ut,C); 
    //testGpuAllocators(&ut,D); 
    const size_t n = 10;
    std::vector<size_t> v1 = {n};
    std::vector<size_t> v2 = {n,n};
    std::vector<size_t> v3 = {n,n,n};
    std::vector<size_t> v4 = {n,n,n,n};
    std::vector<size_t> vLarge = {1e16};
    std::vector<size_t> vEmpty = {0,0,0,0};
    //allocate a simple array
    A.resize(v1); 
    A.resize(v4); 
    A.resize(v3);
    A.resize(v2);
    KernelWrapper<double> K;
    K.setData(A.data(), 4.0, A.length() );
    K.opData(A.data(), A.length() );
    for(size_t i =0; i < A.length(); i++){
        std::cout<< A.data()[i] <<std::endl;
    }   
//make allocator and call this directly, 
//    A.resize(vLarge);
//    A.resize(vEmpty);
    //B.resize(v1);
    AMP::AMPManager::shutdown();
} 
