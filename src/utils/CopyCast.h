#ifndef included_AMP_CopyCast_H_
#define included_AMP_CopyCast_H_

#include <cstddef>

namespace AMP::Utilities {

/*!
 * Helper function to copy and cast (single<->double precision) values between two arrays
 * @param[in]    len      Length of above vectors
 * @param[in]    vec_in   The incoming vector to get the values from
 * @param[inout] vec_out  The outgoing vector to with the up/down-casted values from vec_in
 *                        It is assumed that vec_out is properly allocated
 */
template<typename T1, typename T2>
void copyCast( size_t len, const T1 *vec_in, T2 *vec_out );

} // namespace AMP::Utilities

#endif
