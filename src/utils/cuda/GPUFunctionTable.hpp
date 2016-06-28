#ifndef included_GPUFunctionTable_HPP_
#define included_GPUFunctionTable_HPP_

//Most of this file should be wrapper calls to kernels!!

template <class TYPE, class FUN, class ALLOC>
void GPUFunctionTable::rand(Array<TYPE, FUN, ALLOC> &A)
{
    
