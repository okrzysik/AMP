#include "CopyCast.h"
#include "CopyCast.hpp"

namespace AMP::Utilities {

template void copyCast<double, float>( size_t len, const double *vec_in, float *vec_out );
template void copyCast<float, double>( size_t len, const float *vec_in, double *vec_out );
template void copyCast<float, float>( size_t len, const float *vec_in, float *vec_out );
template void copyCast<double, double>( size_t len, const double *vec_in, double *vec_out );

} // namespace AMP::Utilities