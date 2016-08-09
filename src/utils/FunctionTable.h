#ifndef included_FunctionTable
#define included_FunctionTable


#include "utils/Array.h"

#include <functional>


namespace AMP {


/*!
 * Class FunctionTable is a serial function table class that defines
 *   a series of operations that can be performed on the Array class.
 *   Users can impliment additional versions of the function table that match
 *   the interface to change the behavior of the array class.
 */
class FunctionTable final
{
public:
    /*!
     * Initialize the array with random values
     * @param[in] x         The array to operate on
     */
    template <class TYPE, class FUN>
    static void rand( Array<TYPE, FUN> &x );

    /*!
     * Perform a reduce operator y = f(x)
     * @param[in] op            The function operation
     *                          Note: the operator is a template parameter to improve performance
     * @param[in] x             The array to operate on
     * @param[in] initialValue  The initial value for the reduction (0 for sum, +/- inf for min/max, ...)
     * @return                  The reduction
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline TYPE reduce( LAMBDA &op, const Array<TYPE, FUN> &A, const TYPE& initialValue );

    /*!
     * Perform a reduce operator z = f(x,y)
     * @param[in] op            The function operation
     *                          Note: the operator is a template parameter to improve performance
     * @param[in] x             The first array to operate on
     * @param[in] y             The second array to operate on
     * @param[in] initialValue  The initial value for the reduction (0 for sum, +/- inf for min/max, ...)
     * @return                  The reduction
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline TYPE reduce( LAMBDA &op, const Array<TYPE, FUN> &A, const Array<TYPE, FUN> &B, const TYPE& initialValue );

    /*!
     * Perform a element-wise operation y = f(x)
     * @param[in] fun           The function operation
     *                          Note: the function is a template parameter to improve performance
     * @param[in,out] x         The array to operate on
     * @param[out] y            The output array
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline void transform( LAMBDA &fun, const Array<TYPE, FUN> &x, Array<TYPE, FUN> &y );

    /*!
     * Perform a element-wise operation z = f(x,y)
     * @param[in] fun           The function operation
     *                          Note: the function is a template parameter to improve performance
     * @param[in] x             The first array
     * @param[in] y             The second array
     */
    template <class TYPE, class FUN, typename LAMBDA>
    static inline void transform( LAMBDA &fun,
                                  const Array<TYPE, FUN> &x,
                                  const Array<TYPE, FUN> &y,
                                  Array<TYPE, FUN> &z );

    /*!
     * Multiply two arrays
     * @param[in] a             The first array
     * @param[in] b             The second array
     * @param[out] c            The output array
     */
    template <class TYPE, class FUN>
    static void multiply( const Array<TYPE, FUN> &a, const Array<TYPE, FUN> &b, Array<TYPE, FUN> &c );

    /*!
     * Check if two arrays are approximately equal
     * @param[in] a             The first array
     * @param[in] b             The second array
     * @param[in] TOL           The tolerance
     */
    template <class TYPE, class FUN>
    static bool  equals( const Array<TYPE, FUN> &A, const Array<TYPE, FUN> &B, TYPE tol );

    template <class TYPE>
    static inline void gemmWrapper(char TRANSA, char TRANSB, int M, int N, int K, TYPE alpha, const TYPE* A, int LDA, const TYPE* B, int LDB, TYPE beta, TYPE* C, int LDC);


    /* Specialized Functions */
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

private:
    FunctionTable();

    template <class T>
    static inline void rand( size_t N, T *x );
};


} // namespace AMP

#include "utils/FunctionTable.hpp"

#endif
