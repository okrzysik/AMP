#ifndef included_AMP_Matrix_GetRowHelper_hpp
#define included_AMP_Matrix_GetRowHelper_hpp

#include "AMP/matrices/GetRowHelper.h"

#include <type_traits>


namespace AMP::LinearAlgebra {


/****************************************************************
 * Get rows / size of rows                                       *
 ****************************************************************/
template<class BIGINT_TYPE, class INT_TYPE>
void GetRowHelper::NNZ( BIGINT_TYPE start,
                        BIGINT_TYPE end,
                        INT_TYPE *N_local,
                        INT_TYPE *N_remote ) const
{
    BIGINT_TYPE row = start;
    INT_TYPE i      = 0;
    for ( ; row < end; i++, row++ ) {
        auto rnnz = NNZ( row );
        if ( N_local ) {
            N_local[i] = static_cast<INT_TYPE>( rnnz[0] );
        }
        if ( N_remote ) {
            N_remote[i] = static_cast<INT_TYPE>( rnnz[1] );
        }
    }
}
template<class BIGINT_TYPE, class INT_TYPE>
void GetRowHelper::NNZ( BIGINT_TYPE N_rows,
                        BIGINT_TYPE *rows,
                        INT_TYPE *N_local,
                        INT_TYPE *N_remote ) const
{
    for ( BIGINT_TYPE i = 0; i < N_rows; i++ ) {
        auto rnnz = NNZ( rows[i] );
        if ( N_local ) {
            N_local[i] = static_cast<INT_TYPE>( rnnz[0] );
        }
        if ( N_remote ) {
            N_remote[i] = static_cast<INT_TYPE>( rnnz[1] );
        }
    }
}
template<class BIGINT_TYPE>
void GetRowHelper::getRow( BIGINT_TYPE row, BIGINT_TYPE *local, BIGINT_TYPE *remote ) const
{
    if constexpr ( std::is_integral_v<BIGINT_TYPE> && sizeof( BIGINT_TYPE ) == sizeof( size_t ) ) {
        getRow2( row, reinterpret_cast<size_t *>( local ), reinterpret_cast<size_t *>( remote ) );
    } else {
        auto [N_local, N_remote] = NNZ( row );
        auto local2              = new size_t[N_local];
        auto remote2             = new size_t[N_remote];
        getRow2( row, local2, remote2 );
        if ( local ) {
            for ( size_t i = 0; i < N_local; i++ )
                local[i] = static_cast<BIGINT_TYPE>( local2[i] );
        }
        if ( remote ) {
            for ( size_t i = 0; i < N_remote; i++ )
                remote[i] = static_cast<BIGINT_TYPE>( remote2[i] );
        }
        delete[] local2;
        delete[] remote2;
    }
}
template<class BIGINT_TYPE>
void GetRowHelper::getRow( BIGINT_TYPE start,
                           BIGINT_TYPE end,
                           BIGINT_TYPE **local,
                           BIGINT_TYPE **remote ) const
{
    if constexpr ( std::is_integral_v<BIGINT_TYPE> && sizeof( BIGINT_TYPE ) == sizeof( size_t ) ) {
        const auto end_st = static_cast<size_t>( end );
        for ( size_t row = start, i = 0; row < end_st; row++, i++ )
            getRow2( row,
                     reinterpret_cast<size_t *>( local[i] ),
                     reinterpret_cast<size_t *>( remote[i] ) );
    } else {
        size_t N_max[2] = { 0, 0 };
        for ( size_t row = start, i = 0; row < end; row++, i++ ) {
            auto [N_local, N_remote] = NNZ( row );
            N_max[0]                 = std::max( N_local, N_max[0] );
            N_max[1]                 = std::max( N_remote, N_max[1] );
        }
        auto local2  = new size_t[N_max[0]];
        auto remote2 = new size_t[N_max[1]];
        for ( size_t row = start, i = 0; row < end; row++, i++ ) {
            auto [N_local, N_remote] = NNZ( row );
            getRow2( row, local2, remote2 );
            for ( size_t j = 0; j < N_local; j++ )
                local[i][j] = static_cast<BIGINT_TYPE>( local2[j] );
            for ( size_t j = 0; j < N_remote; j++ )
                remote[i][j] = static_cast<BIGINT_TYPE>( remote2[j] );
        }
        delete[] local2;
        delete[] remote2;
    }
}
template<class BIGINT_TYPE>
void GetRowHelper::getRow( BIGINT_TYPE N_rows,
                           BIGINT_TYPE *rows,
                           BIGINT_TYPE **local,
                           BIGINT_TYPE **remote ) const
{
    if constexpr ( std::is_integral_v<BIGINT_TYPE> && sizeof( BIGINT_TYPE ) == sizeof( size_t ) ) {
        for ( size_t i = i = 0; i < N_rows; i++ )
            getRow2( rows[i],
                     reinterpret_cast<size_t *>( local[i] ),
                     reinterpret_cast<size_t *>( remote[i] ) );
    } else {
        size_t N_max[2] = { 0, 0 };
        for ( size_t i = i = 0; i < N_rows; i++ ) {
            auto [N_local, N_remote] = NNZ( rows[i] );
            N_max[0]                 = std::max( N_local, N_max[0] );
            N_max[1]                 = std::max( N_remote, N_max[1] );
        }
        auto local2  = new size_t[N_max[0]];
        auto remote2 = new size_t[N_max[1]];
        for ( size_t i = i = 0; i < N_rows; i++ ) {
            auto [N_local, N_remote] = NNZ( rows[i] );
            getRow2( rows[i], local2, remote2 );
            for ( size_t j = 0; j < N_local; j++ )
                local[i][j] = static_cast<BIGINT_TYPE>( local2[j] );
            for ( size_t j = 0; j < N_remote; j++ )
                remote[i][j] = static_cast<BIGINT_TYPE>( remote2[j] );
        }
        delete[] local2;
        delete[] remote2;
    }
}


} // namespace AMP::LinearAlgebra

#endif
