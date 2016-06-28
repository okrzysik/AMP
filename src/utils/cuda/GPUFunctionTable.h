#ifndef included_GPUFunctionTable_H_
#define included_GPUFunctionTable_H_

#include "utils/Array.h"


class GPUFunctionTable final
{
public:
    template <class TYPE, class FUN, class ALLOC>
    static void rand(Array<TYPE, FUN, ALLOC> &x);

    template <class TYPE, class FUN, class ALLOC, typename LAMBDA>
    static inline TYPE reduce(LAMBDA &op, const Array<TYPE, FUN, ALLOC>,&A);

    template <class TYPE, class FUN, class Alloc,typename LAMBDA>
    static inline TYPE sum( LAMBDA &op, const Array<TYPE, FUN, ALLOC> &A );

    template <class TYPE, class FUN, class Alloc,typename LAMBDA>
    static inline TYPE reduce( LAMBDA &op, const Array<TYPE, FUN, ALLOC> &A, const Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN,class ALLOC, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN, ALLOC> &x, Array<TYPE, FUN, ALLOC> &y );

    template <class TYPE, class FUN, class ALLOC, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN, ALLOC> &x, const Array<TYPE, FUN, ALLOC> &y, Array<TYPE, FUN, ALLOC> &z );

    template <class TYPE, class FUN>
    static void multiply( const Array<TYPE, FUN, ALLOC> &a, const Array<TYPE, FUN, ALLOC> &b, Array<TYPE, FUN, ALLOC> &c ); 
 
template <class TYPE, class FUN, class ALLOC>
    static bool equals( const Array<TYPE, FUN, ALLOC> &A, const Array<TYPE, FUN, ALLOC> &B, TYPE tol );

private:
    GPUFunctionTable();
  
}; 

#include "utils/FunctionTable.hpp"

#endif
