#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )

    #include "Kokkos_Core.hpp"

namespace AMP::LinearAlgebra {

template<typename Policy, typename Allocator>
auto wrapCSRDiagDataKokkos( CSRMatrixData<Policy, Allocator> *csrData )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const lidx_t nrows                 = static_cast<lidx_t>( csrData->numLocalRows() );
    const lidx_t nnz_tot               = csrData->numberOfNonZerosDiag();
    auto [nnz, cols, cols_loc, coeffs] = csrData->getCSRDiagData();
    auto *rowstarts                    = csrData->getDiagRowStarts();

    // coeffs not marked const so that setScalar and similar will work
    return std::make_tuple(
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( nnz, nrows ),
        Kokkos::View<const gidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( cols, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( cols_loc, nnz_tot ),
        Kokkos::View<scalar_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( coeffs, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( rowstarts,
                                                                                nrows ) );
}

template<typename Policy, typename Allocator>
auto wrapCSROffDiagDataKokkos( CSRMatrixData<Policy, Allocator> *csrData )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const lidx_t nrows                 = static_cast<lidx_t>( csrData->numLocalRows() );
    const lidx_t nnz_tot               = csrData->numberOfNonZerosOffDiag();
    auto [nnz, cols, cols_loc, coeffs] = csrData->getCSROffDiagData();
    auto *rowstarts                    = csrData->getOffDiagRowStarts();

    return std::make_tuple(
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( nnz, nrows ),
        Kokkos::View<const gidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( cols, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( cols_loc, nnz_tot ),
        Kokkos::View<scalar_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( coeffs, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, Kokkos::SharedSpace>( rowstarts,
                                                                                nrows ) );
}

namespace CSRMatOpsKokkosFunctor {

// This functor is based on the one inside KokkosKernels
// Modifications are made to handle our data structures,
// and more importantly our need to handle distributed matrices
template<typename Policy,
         class ExecSpace,
         class IAView,
         class RSView,
         class JAView,
         class AAView,
         class XView,
         class YView>
struct Mult {
    typedef typename Policy::lidx_t lidx_t;
    typedef typename Policy::scalar_t scalar_t;

    lidx_t num_rows;
    lidx_t num_rows_team;
    IAView nnz;
    RSView rowstarts;
    JAView cols_loc;
    AAView coeffs;
    XView inDataBlock;
    YView outDataBlock;

    Mult( lidx_t num_rows_,
          lidx_t num_rows_team_,
          IAView nnz_,
          RSView rowstarts_,
          JAView cols_loc_,
          AAView coeffs_,
          XView inDataBlock_,
          YView outDataBlock_ )
        : num_rows( num_rows_ ),
          num_rows_team( num_rows_team_ ),
          nnz( nnz_ ),
          rowstarts( rowstarts_ ),
          cols_loc( cols_loc_ ),
          coeffs( coeffs_ ),
          inDataBlock( inDataBlock_ ),
          outDataBlock( outDataBlock_ )
    {
    }

    // Calculate product for a single specific row
    KOKKOS_INLINE_FUNCTION
    void operator()( const lidx_t row ) const
    {
        if ( row >= num_rows ) {
            return;
        }
        scalar_t sum  = 0.0;
        const auto nC = nnz( row );
        const auto rs = rowstarts( row );
        for ( lidx_t c = rs; c < rs + nC; ++c ) {
            const auto cl = cols_loc( c );
            sum += coeffs( c ) * inDataBlock( cl );
        }
        outDataBlock( row ) += sum;
    }

    // process a block of rows hierarchically
    KOKKOS_INLINE_FUNCTION
    void operator()( const typename Kokkos::TeamPolicy<ExecSpace>::member_type &tm ) const
    {
        const auto lRank = tm.league_rank();
        const auto fRow  = lRank * num_rows_team;
        Kokkos::parallel_for( Kokkos::TeamThreadRange( tm, num_rows_team ),
                              [&]( const lidx_t tIdx ) {
                                  const auto row = fRow + tIdx;
                                  if ( row >= num_rows ) {
                                      return;
                                  }
                                  scalar_t sum  = 0.0;
                                  const auto nC = nnz( row );
                                  const auto rs = rowstarts( row );
                                  Kokkos::parallel_reduce(
                                      Kokkos::ThreadVectorRange( tm, nC ),
                                      [&]( lidx_t &c, scalar_t &lsum ) {
                                          const auto cl = cols_loc( rs + c );
                                          lsum += coeffs( rs + c ) * inDataBlock( cl );
                                      },
                                      sum );
                                  outDataBlock( row ) += sum;
                              } );
    }
};

template<typename Policy,
         class ExecSpace,
         class IAView,
         class RSView,
         class JAView,
         class AAView,
         class XView,
         class YView>
struct MultTranspose {
    typedef typename Policy::lidx_t lidx_t;
    typedef typename Policy::scalar_t scalar_t;

    lidx_t num_rows;
    lidx_t num_rows_team;
    IAView nnz;
    RSView rowstarts;
    JAView cols_loc;
    AAView coeffs;
    XView inDataBlock;
    YView outDataBlock;

    MultTranspose( lidx_t num_rows_,
                   lidx_t num_rows_team_,
                   IAView nnz_,
                   RSView rowstarts_,
                   JAView cols_loc_,
                   AAView coeffs_,
                   XView inDataBlock_,
                   YView outDataBlock_ )
        : num_rows( num_rows_ ),
          num_rows_team( num_rows_team_ ),
          nnz( nnz_ ),
          rowstarts( rowstarts_ ),
          cols_loc( cols_loc_ ),
          coeffs( coeffs_ ),
          inDataBlock( inDataBlock_ ),
          outDataBlock( outDataBlock_ )
    {
        // TODO: Add assertion that YView has atomic memory trait
    }

    // Calculate product for a single specific row
    KOKKOS_INLINE_FUNCTION
    void operator()( const lidx_t row ) const
    {
        if ( row >= num_rows ) {
            return;
        }
        const auto nC = nnz( row );
        const auto rs = rowstarts( row );
        const auto xi = inDataBlock( row );
        for ( lidx_t c = rs; c < rs + nC; ++c ) {
            const auto cl = cols_loc( c );
            outDataBlock( cl ) += xi * coeffs( c );
        }
    }

    // process a block of rows hierarchically
    KOKKOS_INLINE_FUNCTION
    void operator()( const typename Kokkos::TeamPolicy<ExecSpace>::member_type &tm ) const
    {
        const auto lRank = tm.league_rank();
        const auto fRow  = lRank * num_rows_team;
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange( tm, num_rows_team ), [&]( const lidx_t tIdx ) {
                const auto row = fRow + tIdx;
                if ( row >= num_rows ) {
                    return;
                }
                const auto nC = nnz( row );
                const auto rs = rowstarts( row );
                const auto xi = inDataBlock( row );
                Kokkos::parallel_for( Kokkos::ThreadVectorRange( tm, nC ), [&]( lidx_t &c ) {
                    const auto cl = cols_loc( rs + c );
                    outDataBlock( cl ) += xi * coeffs( rs + c );
                } );
            } );
    }
};

} // namespace CSRMatOpsKokkosFunctor

} // namespace AMP::LinearAlgebra
