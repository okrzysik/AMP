#ifndef included_AMP_CopyCast_kokkos_HPP_
#define included_AMP_CopyCast_kokkos_HPP_

#include "AMP/AMP_TPLs.h"
#include <Kokkos_Core.hpp>
#include <Kokkos_Macros.hpp>

#include "AMP/utils/memory.h"
#include <memory>

namespace AMP::Utilities {
#ifdef USE_OPENMP
typedef Kokkos::OpenMP HostExecSpace;
#else
typedef Kokkos::Serial HostExecSpace;
#endif
typedef Kokkos::RangePolicy<HostExecSpace> host_range_policy;

/*!
 * Helper function to copy and cast (single<->double precision) values between two arrays
 * @param[in]    len      Length of above vectors
 * @param[in]    vec_in   The incoming vector to get the values from
 * @param[inout] vec_out  The outgoing vector to with the up/down-casted values from vec_in
 *                        It is assumed that vec_out is properly allocated
 */
template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::host> {
    void operator()( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        Kokkos::parallel_for(
            "Copy cast", host_range_policy( 0, len ), KOKKOS_LAMBDA( const size_t i ) {
                AMP_ASSERT( std::abs( vec_in[i] ) <= std::numeric_limits<T2>::max() );
                vec_out[i] = static_cast<T2>( vec_in[i] );
            } );
    }
};

template<typename T1, typename T2>
struct copyCast_<T1, T2, AMP::Utilities::MemoryType::unregistered> {
    void operator()( size_t len, const T1 *vec_in, T2 *vec_out )
    {
        Kokkos::parallel_for(
            "Copy cast", host_range_policy( 0, len ), KOKKOS_LAMBDA( const size_t i ) {
                AMP_ASSERT( std::abs( vec_in[i] ) <= std::numeric_limits<T2>::max() );
                vec_out[i] = static_cast<T2>( vec_in[i] );
            } );
    }
};

} // namespace AMP::Utilities

// #include "CopyCast.hpp"

#endif
