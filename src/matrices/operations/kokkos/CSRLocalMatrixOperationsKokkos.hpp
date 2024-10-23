#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRLocalMatrixData.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/kokkos/CSRLocalMatrixOperationsKokkos.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/memory.h"
#include "AMP/vectors/Vector.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )

    #include "Kokkos_Core.hpp"

namespace AMP::LinearAlgebra {

template<typename Policy, typename Allocator, class ViewSpace>
auto wrapCSRDataKokkos( std::shared_ptr<CSRLocalMatrixData<Policy, Allocator>> A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const lidx_t nrows                 = static_cast<lidx_t>( A->numLocalRows() );
    const lidx_t nnz_tot               = A->numberOfNonZeros();
    auto [nnz, cols, cols_loc, coeffs] = A->getDataFields();
    auto *rowstarts                    = A->getRowStarts();

    // coeffs not marked const so that setScalar and similar will work
    return std::make_tuple(
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( nnz, nrows ),
        Kokkos::View<const gidx_t *, Kokkos::LayoutRight, ViewSpace>( cols, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( cols_loc, nnz_tot ),
        Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace>( coeffs, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( rowstarts, nrows + 1 ) );
}

namespace CSRMatOpsKokkosFunctor {

// This functor is based on the one inside KokkosKernels
// Modifications are made to handle our data structures,
// and more importantly our need to handle distributed matrices
template<typename Policy,
         class ExecSpace,
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
    RSView rowstarts;
    JAView cols_loc;
    AAView coeffs;
    XView inDataBlock;
    YView outDataBlock;

    Mult( lidx_t num_rows_,
          lidx_t num_rows_team_,
          RSView rowstarts_,
          JAView cols_loc_,
          AAView coeffs_,
          XView inDataBlock_,
          YView outDataBlock_ )
        : num_rows( num_rows_ ),
          num_rows_team( num_rows_team_ ),
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
        const auto rs = rowstarts( row );
        const auto re = rowstarts( row + 1 );
        for ( lidx_t c = rs; c < re; ++c ) {
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
                                  const auto rs = rowstarts( row );
                                  const auto re = rowstarts( row + 1 );
                                  Kokkos::parallel_reduce(
                                      Kokkos::ThreadVectorRange( tm, re - rs ),
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
    RSView rowstarts;
    JAView cols_loc;
    AAView coeffs;
    XView inDataBlock;
    YView outDataBlock;

    MultTranspose( lidx_t num_rows_,
                   lidx_t num_rows_team_,
                   RSView rowstarts_,
                   JAView cols_loc_,
                   AAView coeffs_,
                   XView inDataBlock_,
                   YView outDataBlock_ )
        : num_rows( num_rows_ ),
          num_rows_team( num_rows_team_ ),
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
        const auto rs = rowstarts( row );
        const auto re = rowstarts( row + 1 );
        const auto xi = inDataBlock( row );
        for ( lidx_t c = rs; c < re; ++c ) {
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
                const auto rs = rowstarts( row );
                const auto re = rowstarts( row + 1 );
                const auto xi = inDataBlock( row );
                Kokkos::parallel_for( Kokkos::ThreadVectorRange( tm, re - rs ), [&]( lidx_t &c ) {
                    const auto cl = cols_loc( rs + c );
                    outDataBlock( cl ) += xi * coeffs( rs + c );
                } );
            } );
    }
};

} // namespace CSRMatOpsKokkosFunctor

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::mult(
    const typename Policy::scalar_t *in,
    std::shared_ptr<LocalMatrixData> A,
    typename Policy::scalar_t *out )
{
    using scalar_t = typename Policy::scalar_t;
    using lidx_t   = typename Policy::lidx_t;

    const auto nRows = A->numLocalRows();
    const auto nCols = A->numUniqueColumns();

    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 ViewSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        inView( in, nCols );
    Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace> outView( out, nRows );

    const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
        wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );

    const lidx_t team_rows     = 64;
    const lidx_t vector_length = 8;
    const lidx_t num_teams     = ( nRows + team_rows - 1 ) / team_rows;

    CSRMatOpsKokkosFunctor::Mult<Policy,
                                 ExecSpace,
                                 decltype( rowstarts ),
                                 decltype( cols_loc ),
                                 decltype( coeffs ),
                                 decltype( inView ),
                                 decltype( outView )>
        ftor( nRows, team_rows, rowstarts, cols_loc, coeffs, inView, outView );

    if constexpr ( std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace> ) {
        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> team_policy(
            d_exec_space, num_teams, Kokkos::AUTO, vector_length );
        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (local - team)", team_policy, ftor );
    } else {
        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (local - flat)",
                              Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
                              ftor );
    }
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    multTranspose( const typename Policy::scalar_t *in,
                   std::shared_ptr<LocalMatrixData> A,
                   typename Policy::scalar_t *out )
{
    using scalar_t = typename Policy::scalar_t;
    using lidx_t   = typename Policy::lidx_t;

    const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
        wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );

    const auto nRows    = A->numLocalRows();
    const auto nCols    = A->numLocalColumns();
    const auto nColsUnq = A->numUniqueColumns();

    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 ViewSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        inView( in, nCols );

    // vvals only on host, use temporary and copy as needed after
    Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>
        outView( out, nColsUnq );

    const lidx_t team_rows     = 64;
    const lidx_t vector_length = 8;
    const lidx_t num_teams     = ( nRows + team_rows - 1 ) / team_rows;

    CSRMatOpsKokkosFunctor::MultTranspose<Policy,
                                          ExecSpace,
                                          decltype( rowstarts ),
                                          decltype( cols_loc ),
                                          decltype( coeffs ),
                                          decltype( inView ),
                                          decltype( outView )>
        ftor( nRows, team_rows, rowstarts, cols_loc, coeffs, inView, outView );

    if constexpr ( std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace> && false ) {
        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> team_policy(
            d_exec_space, num_teams, Kokkos::AUTO, vector_length );
        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::multTranspose (local - team)", team_policy, ftor );
    } else {
        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::multTranspose (local - flat)",
                              Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
                              ftor );
    }
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    scale( typename Policy::scalar_t alpha, std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t = typename Policy::lidx_t;

    // const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    auto coeffs     = std::get<3>( vTpl );

    const auto tnnz = A->numberOfNonZeros();

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::scale",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz ),
        KOKKOS_LAMBDA( lidx_t n ) { coeffs( n ) *= alpha; } );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    matMultiply( std::shared_ptr<LocalMatrixData>,
                 std::shared_ptr<LocalMatrixData>,
                 std::shared_ptr<LocalMatrixData> )
{
    AMP_WARNING( "SpGEMM for CSRLocalMatrixOperationsKokkos not implemented" );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::axpy(
    typename Policy::scalar_t alpha,
    std::shared_ptr<LocalMatrixData> X,
    std::shared_ptr<LocalMatrixData> Y )
{
    using gidx_t = typename Policy::gidx_t;

    // const auto [nnzX, colsX, cols_locX, coeffsX, rowstartsX] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( X );
    const auto vTplX = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( X );
    auto coeffsX     = std::get<3>( vTplX );

    // const auto [nnzY, colsY, cols_locY, coeffsY, rowstartsY] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( Y );
    const auto vTplY = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( Y );
    auto coeffsY     = std::get<3>( vTplY );

    const auto tnnz = X->numberOfNonZeros();

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::axpy",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz ),
        KOKKOS_LAMBDA( gidx_t n ) { coeffsY( n ) += alpha * coeffsX( n ); } );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    setScalar( typename Policy::scalar_t alpha, std::shared_ptr<LocalMatrixData> A )
{
    const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
        wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    // const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );

    Kokkos::deep_copy( d_exec_space, coeffs, alpha );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::zero(
    std::shared_ptr<LocalMatrixData> A )
{
    setScalar( 0.0, A );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    setDiagonal( const typename Policy::scalar_t *in, std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    const auto beginRow = A->beginRow();
    const auto nRows    = A->numLocalRows();

    // const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    auto nnz        = std::get<0>( vTpl );
    auto cols       = std::get<1>( vTpl );
    auto coeffs     = std::get<3>( vTpl );
    auto rowstarts  = std::get<4>( vTpl );

    Kokkos::View<const scalar_t *, Kokkos::LayoutRight> vvals( in, nRows );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::setDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz( row );
            const auto rs = rowstarts( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols( rs + c ) == static_cast<gidx_t>( row ) + beginRow ) {
                    coeffs( rs + c ) = vvals( row );
                    break;
                }
            }
        } );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    setIdentity( std::shared_ptr<LocalMatrixData> A )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    const auto beginRow = A->beginRow();
    const auto nRows    = A->numLocalRows();

    // const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    auto nnz        = std::get<0>( vTpl );
    auto cols       = std::get<1>( vTpl );
    auto coeffs     = std::get<3>( vTpl );
    auto rowstarts  = std::get<4>( vTpl );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::setDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz( row );
            const auto rs = rowstarts( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols( rs + c ) == static_cast<gidx_t>( row ) + beginRow ) {
                    coeffs( rs + c ) = 1.0;
                    break;
                }
            }
        } );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    extractDiagonal( std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *buf )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;

    const auto beginRow = A->beginRow();
    const auto nRows    = static_cast<lidx_t>( A->numLocalRows() );

    // const auto [nnz, cols, cols_loc, coeffs, rowstarts] =
    //     wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    auto nnz        = std::get<0>( vTpl );
    auto cols       = std::get<1>( vTpl );
    auto coeffs     = std::get<3>( vTpl );
    auto rowstarts  = std::get<4>( vTpl );

    Kokkos::View<scalar_t *, Kokkos::LayoutRight> vvals( buf, nRows );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::extractDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz( row );
            const auto rs = rowstarts( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols( rs + c ) == static_cast<gidx_t>( row ) + beginRow ) {
                    vvals( row ) = coeffs( rs + c );
                    break;
                }
            }
        } );
}

template<typename Policy, class Allocator, class ExecSpace, class ViewSpace, class LocalMatrixData>
void CSRLocalMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace, LocalMatrixData>::
    LinfNorm( std::shared_ptr<LocalMatrixData> A, typename Policy::scalar_t *rowSums ) const
{
    using scalar_t = typename Policy::scalar_t;
    using lidx_t   = typename Policy::lidx_t;

    const auto nRows = A->numLocalRows();

    const auto vTpl = wrapCSRDataKokkos<Policy, Allocator, ViewSpace>( A );
    auto nnz        = std::get<0>( vTpl );
    auto coeffs     = std::get<3>( vTpl );
    auto rowstarts  = std::get<4>( vTpl );

    Kokkos::View<scalar_t *, Kokkos::LayoutRight> sums( rowSums, nRows );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::extractDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz( row );
            const auto rs = rowstarts( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                sums( row ) += Kokkos::fabs( coeffs( rs + c ) );
            }
        } );
}

} // namespace AMP::LinearAlgebra

#endif
