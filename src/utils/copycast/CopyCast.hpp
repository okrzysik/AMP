#ifndef included_AMP_CopyCast_HPP_
#define included_AMP_CopyCast_HPP_

#include "AMP/AMP_TPLs.h"
#include "AMP/utils/memory.h"
#include "AMP/utils/Utilities.h"


#if defined( USE_DEVICE )
    #include "AMP/utils/copycast/device/DevCopyCast.h"
    #include "AMP/utils/copycast/device/DevCopyCast.hpp"
#elif defined( AMP_USE_KOKKOS )
    #include "AMP/utils/copycast/kokkos/CopyCast.h"
    #include "AMP/utils/copycast/kokkos/CopyCast.hpp"
#elif defined( USE_OPENMP )
    #include "AMP/utils/copycast/openmp/CopyCast.h"
    #include "AMP/utils/copycast/openmp/CopyCast.hpp"
#else
    #include "AMP/utils/copycast/serial/CopyCast.h"
    #include "AMP/utils/copycast/serial/CopyCast.hpp"
#endif

#include <iostream>
#include <memory>

namespace AMP::Utilities {

/*!
 * Helper function to copy and cast (single<->double precision) values between two arrays
 * @param[in]    len      Length of above vectors
 * @param[in]    vec_in   The incoming vector to get the values from
 * @param[inout] vec_out  The outgoing vector to with the up/down-casted values from vec_in
 *                        It is assumed that vec_out is properly allocated
 */
template<typename T1, typename T2>
void copyCast( size_t len, const T1 *vec_in, T2 *vec_out )
{
    auto memtype_in  = getMemoryType( vec_in );
    auto memtype_out = getMemoryType( vec_out );
    AMP_ASSERT( memtype_in == memtype_out );
    if ( memtype_in == AMP::Utilities::MemoryType::host ) {
        copyCast_<T1, T2, AMP::Utilities::MemoryType::host> apply;
        apply( len, vec_in, vec_out );
#ifdef USE_DEVICE
    } else if ( memtype_in == AMP::Utilities::MemoryType::managed ) {
        copyCast_<T1, T2, AMP::Utilities::MemoryType::managed> apply;
        apply( len, vec_in, vec_out );
    } else if ( memtype_in == AMP::Utilities::MemoryType::device ) {
        copyCast_<T1, T2, AMP::Utilities::MemoryType::device> apply;
        apply( len, vec_in, vec_out );
#endif
    } else { // unregistered
        copyCast_<T1, T2, AMP::Utilities::MemoryType::unregistered> apply;
        apply( len, vec_in, vec_out );
    }
}

} // namespace AMP::Utilities

#endif
