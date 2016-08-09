#ifndef included_GPUFunctionTable_H_
#define included_GPUFunctionTable_H_

#include "utils/Array.h"

namespace AMP{

class GPUFunctionTable final
{
public:
    template <class TYPE, class FUN, class ALLOC>
    static inline void rand(Array<TYPE, FUN, ALLOC> &x);

    template <class TYPE, class FUN, class ALLOC, typename LAMBDA>
    static inline TYPE reduce(LAMBDA &op, const Array<TYPE, FUN, ALLOC> &A);
    
    template <class TYPE, class FUN, class ALLOC,typename LAMBDA>
    static inline TYPE reduce( LAMBDA &op, const Array<TYPE, FUN, ALLOC> &A, const Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN,class ALLOC, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN, ALLOC> &x, Array<TYPE, FUN, ALLOC> &y );

    template <class TYPE, class FUN, class ALLOC, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN, ALLOC> &x, const Array<TYPE, FUN, ALLOC> &y, Array<TYPE, FUN, ALLOC> &z );

    template <class TYPE, class FUN, class ALLOC>
    static void multiply( const Array<TYPE, FUN, ALLOC> &a, const Array<TYPE, FUN, ALLOC> &b, Array<TYPE, FUN, ALLOC> &c ); 
 
template <class TYPE, class FUN, class ALLOC>
    static bool equals( const Array<TYPE, FUN, ALLOC> &A, const Array<TYPE, FUN, ALLOC> &B, TYPE tol );

    template <class TYPE, class FUN, class ALLOC>
    static void transformReLU(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);
     
    template <class TYPE, class FUN, class ALLOC>
    static void transformAbs(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN, class ALLOC>
    static void transformTanh(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN, class ALLOC>
    static void transformHardTanh(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);
    
    template <class TYPE, class FUN, class ALLOC>
    static void transformSigmoid(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN, class ALLOC>
    static void transformSoftPlus(const Array<TYPE, FUN, ALLOC> &A, Array<TYPE, FUN, ALLOC> &B);

    template <class TYPE, class FUN, class ALLOC>
    static TYPE sum(const Array<TYPE, FUN, ALLOC> &A); 

     
    template <class TYPE>
    static inline void gemmWrapper(char TRANSA, char TRANSB, int M, int N, int K, TYPE alpha, const TYPE* A, int LDA, const TYPE* B, int LDB, TYPE beta, TYPE* C, int LDC);

    GPUFunctionTable(){};
private:

    template <class TYPE>
    static inline void rand(size_t N, TYPE *x);
}; 

}
#include "utils/cuda/GPUFunctionTable.hpp"


#endif
