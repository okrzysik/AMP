#ifndef included_AMP_DelaunayHelpers
#define included_AMP_DelaunayHelpers

#include <stdlib.h>


namespace AMP {


/********************************************************************
 * Solver helper class                                               *
 * Note: The matrix is stored in column-major order                  *
 * Note: For the determinant to be exact, we must support a signed   *
 *   range of +- N^D, where N is the largest input value and D is    *
 *   the size of the matrix.                                         *
 * Note: The solve functions do not normalize the inverse:           *
 *      x = det * M^-1 * b                                           *
 *   Instead they return the normalization constant det so the user  *
 *   can perform the normalization if necessary.  This helps to      *
 *   preserve accuracy.                                              *
 ********************************************************************/
template<int NDIM>
class DelaunayHelpers
{
public:
    // Compute the determinant of the matrix
    template<class TYPE>
    static inline TYPE det( const TYPE *M );
    // Solve the linear system A*x = d*b, where d is the determinant of the matrix
    // Note:  This can be done exactly, while solving A*x = b can not
    template<class TYPE>
    static inline void solve_system( const TYPE *M, const TYPE *b, TYPE *x, TYPE &det );
};
template<int NDIM>
template<class TYPE>
inline TYPE DelaunayHelpers<NDIM>::det( const TYPE *M )
{
    if constexpr ( NDIM == 1 ) {
        return M[0];
    } else if constexpr ( NDIM == 2 ) {
        return M[0] * M[3] - M[1] * M[2];
    } else if constexpr ( NDIM == 3 ) {
        TYPE det( 0 );
        det += M[0] * ( M[4] * M[8] - M[7] * M[5] );
        det -= M[3] * ( M[1] * M[8] - M[7] * M[2] );
        det += M[6] * ( M[1] * M[5] - M[4] * M[2] );
        return det;
    } else if constexpr ( NDIM == 4 ) {
        TYPE tmp[6];
        tmp[0] = M[2] * M[7] - M[6] * M[3];
        tmp[1] = M[2] * M[11] - M[10] * M[3];
        tmp[2] = M[2] * M[15] - M[14] * M[3];
        tmp[3] = M[6] * M[11] - M[10] * M[7];
        tmp[4] = M[6] * M[15] - M[14] * M[7];
        tmp[5] = M[10] * M[15] - M[14] * M[11];
        TYPE det( 0 );
        det += M[0] * ( M[5] * tmp[5] - M[9] * tmp[4] + M[13] * tmp[3] );
        det -= M[4] * ( M[1] * tmp[5] - M[9] * tmp[2] + M[13] * tmp[1] );
        det += M[8] * ( M[1] * tmp[4] - M[5] * tmp[2] + M[13] * tmp[0] );
        det -= M[12] * ( M[1] * tmp[3] - M[5] * tmp[1] + M[9] * tmp[0] );
        return det;
    } else {
        throw std::logic_error( "Not programmed" );
    }
}
template<int NDIM>
template<class TYPE>
inline void
DelaunayHelpers<NDIM>::solve_system( const TYPE *M, const TYPE *b, TYPE *x, TYPE &det_M )
{
    if constexpr ( NDIM == 1 ) {
        det_M = M[0];
        x[0]  = b[0];
    } else if constexpr ( NDIM == 2 ) {
        det_M = M[0] * M[3] - M[1] * M[2];
        x[0]  = M[3] * b[0] - M[2] * b[1];
        x[1]  = M[0] * b[1] - M[1] * b[0];
    } else if constexpr ( NDIM == 3 ) {
        TYPE inv[9];
        inv[0] = M[4] * M[8] - M[7] * M[5];
        inv[1] = M[7] * M[2] - M[1] * M[8];
        inv[2] = M[1] * M[5] - M[4] * M[2];
        inv[3] = M[6] * M[5] - M[3] * M[8];
        inv[4] = M[0] * M[8] - M[6] * M[2];
        inv[5] = M[3] * M[2] - M[0] * M[5];
        inv[6] = M[3] * M[7] - M[6] * M[4];
        inv[7] = M[6] * M[1] - M[0] * M[7];
        inv[8] = M[0] * M[4] - M[3] * M[1];
        det_M  = M[0] * inv[0] + M[3] * inv[1] + M[6] * inv[2];
        x[0]   = inv[0] * b[0] + inv[3] * b[1] + inv[6] * b[2];
        x[1]   = inv[1] * b[0] + inv[4] * b[1] + inv[7] * b[2];
        x[2]   = inv[2] * b[0] + inv[5] * b[1] + inv[8] * b[2];
    } else if constexpr ( NDIM == 4 ) {
        TYPE inv[16];
        TYPE m00 = M[0], m01 = M[4], m02 = M[8], m03 = M[12];
        TYPE m10 = M[1], m11 = M[5], m12 = M[9], m13 = M[13];
        TYPE m20 = M[2], m21 = M[6], m22 = M[10], m23 = M[14];
        TYPE m30 = M[3], m31 = M[7], m32 = M[11], m33 = M[15];
        det_M  = det( M );
        inv[0] = m12 * m23 * m31 - m13 * m22 * m31 + m13 * m21 * m32 - m11 * m23 * m32 -
                 m12 * m21 * m33 + m11 * m22 * m33;
        inv[1] = m03 * m22 * m31 - m02 * m23 * m31 - m03 * m21 * m32 + m01 * m23 * m32 +
                 m02 * m21 * m33 - m01 * m22 * m33;
        inv[2] = m02 * m13 * m31 - m03 * m12 * m31 + m03 * m11 * m32 - m01 * m13 * m32 -
                 m02 * m11 * m33 + m01 * m12 * m33;
        inv[3] = m03 * m12 * m21 - m02 * m13 * m21 - m03 * m11 * m22 + m01 * m13 * m22 +
                 m02 * m11 * m23 - m01 * m12 * m23;
        inv[4] = m13 * m22 * m30 - m12 * m23 * m30 - m13 * m20 * m32 + m10 * m23 * m32 +
                 m12 * m20 * m33 - m10 * m22 * m33;
        inv[5] = m02 * m23 * m30 - m03 * m22 * m30 + m03 * m20 * m32 - m00 * m23 * m32 -
                 m02 * m20 * m33 + m00 * m22 * m33;
        inv[6] = m03 * m12 * m30 - m02 * m13 * m30 - m03 * m10 * m32 + m00 * m13 * m32 +
                 m02 * m10 * m33 - m00 * m12 * m33;
        inv[7] = m02 * m13 * m20 - m03 * m12 * m20 + m03 * m10 * m22 - m00 * m13 * m22 -
                 m02 * m10 * m23 + m00 * m12 * m23;
        inv[8] = m11 * m23 * m30 - m13 * m21 * m30 + m13 * m20 * m31 - m10 * m23 * m31 -
                 m11 * m20 * m33 + m10 * m21 * m33;
        inv[9] = m03 * m21 * m30 - m01 * m23 * m30 - m03 * m20 * m31 + m00 * m23 * m31 +
                 m01 * m20 * m33 - m00 * m21 * m33;
        inv[10] = m01 * m13 * m30 - m03 * m11 * m30 + m03 * m10 * m31 - m00 * m13 * m31 -
                  m01 * m10 * m33 + m00 * m11 * m33;
        inv[11] = m03 * m11 * m20 - m01 * m13 * m20 - m03 * m10 * m21 + m00 * m13 * m21 +
                  m01 * m10 * m23 - m00 * m11 * m23;
        inv[12] = m12 * m21 * m30 - m11 * m22 * m30 - m12 * m20 * m31 + m10 * m22 * m31 +
                  m11 * m20 * m32 - m10 * m21 * m32;
        inv[13] = m01 * m22 * m30 - m02 * m21 * m30 + m02 * m20 * m31 - m00 * m22 * m31 -
                  m01 * m20 * m32 + m00 * m21 * m32;
        inv[14] = m02 * m11 * m30 - m01 * m12 * m30 - m02 * m10 * m31 + m00 * m12 * m31 +
                  m01 * m10 * m32 - m00 * m11 * m32;
        inv[15] = m01 * m12 * m20 - m02 * m11 * m20 + m02 * m10 * m21 - m00 * m12 * m21 -
                  m01 * m10 * m22 + m00 * m11 * m22;
        x[0] = ( inv[0] * b[0] + inv[1] * b[1] + inv[2] * b[2] + inv[3] * b[3] ) / det_M;
        x[1] = ( inv[4] * b[0] + inv[5] * b[1] + inv[6] * b[2] + inv[7] * b[3] ) / det_M;
        x[2] = ( inv[8] * b[0] + inv[9] * b[1] + inv[10] * b[2] + inv[11] * b[3] ) / det_M;
        x[3] = ( inv[12] * b[0] + inv[13] * b[1] + inv[14] * b[2] + inv[15] * b[3] ) / det_M;
    }
}


} // namespace AMP

#endif
