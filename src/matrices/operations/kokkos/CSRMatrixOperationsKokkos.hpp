#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/kokkos/CSRMatrixOperationsKokkos.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )

    #include "Kokkos_Core.hpp"

namespace AMP::LinearAlgebra {

template<typename Policy, typename Allocator, class ViewSpace>
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
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( nnz, nrows ),
        Kokkos::View<const gidx_t *, Kokkos::LayoutRight, ViewSpace>( cols, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( cols_loc, nnz_tot ),
        Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace>( coeffs, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( rowstarts, nrows ) );
}

template<typename Policy, typename Allocator, class ViewSpace>
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
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( nnz, nrows ),
        Kokkos::View<const gidx_t *, Kokkos::LayoutRight, ViewSpace>( cols, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( cols_loc, nnz_tot ),
        Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace>( coeffs, nnz_tot ),
        Kokkos::View<const lidx_t *, Kokkos::LayoutRight, ViewSpace>( rowstarts, nrows ) );
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

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::mult(
    std::shared_ptr<const Vector> in, MatrixData const &A, std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsKokkos::mult" );
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    const auto nCols = static_cast<lidx_t>( csrData->numLocalColumns() );

    auto inData        = in->getVectorData();
    const auto &ghosts = inData->getGhosts();
    const auto nGhosts = static_cast<lidx_t>( ghosts.size() );
    auto outData       = out->getVectorData();

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    AMP_DEBUG_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with one data block" );

    AMP_DEBUG_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with accessible ghosts" );


    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 ViewSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        inDataBlock( inData->getRawDataBlock<scalar_t>( 0 ), nCols );
    Kokkos::View<scalar_t *, Kokkos::LayoutRight, ViewSpace> outDataBlock(
        outData->getRawDataBlock<scalar_t>( 0 ), nRows );
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 ViewSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        ghostDataBlock( ghosts.data(), nGhosts );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl        = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_d       = std::get<0>( vtpl );
        auto cols_loc_d  = std::get<2>( vtpl );
        auto coeffs_d    = std::get<3>( vtpl );
        auto rowstarts_d = std::get<4>( vtpl );

        const lidx_t team_rows     = 64;
        const lidx_t vector_length = 8;
        const lidx_t num_teams     = ( nRows + team_rows - 1 ) / team_rows;

        CSRMatOpsKokkosFunctor::Mult<Policy,
                                     ExecSpace,
                                     decltype( nnz_d ),
                                     decltype( rowstarts_d ),
                                     decltype( cols_loc_d ),
                                     decltype( coeffs_d ),
                                     decltype( inDataBlock ),
                                     decltype( outDataBlock )>
            ftor( nRows,
                  team_rows,
                  nnz_d,
                  rowstarts_d,
                  cols_loc_d,
                  coeffs_d,
                  inDataBlock,
                  outDataBlock );

        if constexpr ( std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace> ) {
            Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> team_policy(
                d_exec_space, num_teams, Kokkos::AUTO, vector_length );
            Kokkos::parallel_for(
                "CSRMatrixOperationsKokkos::mult (local - team)", team_policy, ftor );
        } else {
            Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (local - flat)",
                                  Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
                                  ftor );
        }
    }

    if ( csrData->hasOffDiag() ) {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl         = wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_od       = std::get<0>( vtpl );
        auto cols_loc_od  = std::get<2>( vtpl );
        auto coeffs_od    = std::get<3>( vtpl );
        auto rowstarts_od = std::get<4>( vtpl );

        CSRMatOpsKokkosFunctor::Mult<Policy,
                                     ExecSpace,
                                     decltype( nnz_od ),
                                     decltype( rowstarts_od ),
                                     decltype( cols_loc_od ),
                                     decltype( coeffs_od ),
                                     decltype( ghostDataBlock ),
                                     decltype( outDataBlock )>
            ftor( nRows,
                  1,
                  nnz_od,
                  rowstarts_od,
                  cols_loc_od,
                  coeffs_od,
                  ghostDataBlock,
                  outDataBlock );

        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (ghost)",
                              Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
                              ftor );
    }

    d_exec_space.fence(); // get rid of this eventually
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::multTranspose(
    std::shared_ptr<const Vector> in, MatrixData const &A, std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsKokkos::multTranspose" );

    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    const auto nCols = static_cast<lidx_t>( csrData->numLocalColumns() );
    auto inData      = in->getVectorData();
    auto outData     = out->getVectorData();

    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *, Kokkos::LayoutRight, ViewSpace> inDataBlock(
        inData->getRawDataBlock<scalar_t>( 0 ), nRows );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl        = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_d       = std::get<0>( vtpl );
        auto cols_loc_d  = std::get<2>( vtpl );
        auto coeffs_d    = std::get<3>( vtpl );
        auto rowstarts_d = std::get<4>( vtpl );

        // Make temporary views for output columns and values
        auto outDataBlock = Kokkos::
            View<scalar_t *, Kokkos::LayoutRight, ViewSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>(
                outData->getRawDataBlock<scalar_t>( 0 ), nCols );

        const lidx_t team_rows     = 64;
        const lidx_t vector_length = 8;
        const lidx_t num_teams     = ( nRows + team_rows - 1 ) / team_rows;

        CSRMatOpsKokkosFunctor::MultTranspose<Policy,
                                              ExecSpace,
                                              decltype( nnz_d ),
                                              decltype( rowstarts_d ),
                                              decltype( cols_loc_d ),
                                              decltype( coeffs_d ),
                                              decltype( inDataBlock ),
                                              decltype( outDataBlock )>
            ftor( nRows,
                  team_rows,
                  nnz_d,
                  rowstarts_d,
                  cols_loc_d,
                  coeffs_d,
                  inDataBlock,
                  outDataBlock );

        if constexpr ( std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace> ) {
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

    if ( csrData->hasOffDiag() ) {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl         = wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_od       = std::get<0>( vtpl );
        auto cols_loc_od  = std::get<2>( vtpl );
        auto coeffs_od    = std::get<3>( vtpl );
        auto rowstarts_od = std::get<4>( vtpl );

        // get off diag map
        auto colMap        = csrData->getOffDiagColumnMap();
        const auto num_unq = csrData->getOffDiagColumnMapSize();

        std::vector<size_t> rcols( num_unq );
        std::transform(
            colMap, colMap + num_unq, rcols.begin(), []( size_t col ) -> size_t { return col; } );

        // Make temporary view for output values
        auto vvals = Kokkos::
            View<scalar_t *, Kokkos::LayoutRight, ViewSpace, Kokkos::MemoryTraits<Kokkos::Atomic>>(
                "vvals", num_unq );
        Kokkos::deep_copy( d_exec_space, vvals, 0.0 );

        CSRMatOpsKokkosFunctor::MultTranspose<Policy,
                                              ExecSpace,
                                              decltype( nnz_od ),
                                              decltype( rowstarts_od ),
                                              decltype( cols_loc_od ),
                                              decltype( coeffs_od ),
                                              decltype( inDataBlock ),
                                              decltype( vvals )>
            ftor( nRows, 1, nnz_od, rowstarts_od, cols_loc_od, coeffs_od, inDataBlock, vvals );

        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::multTranspose (ghost)",
                              Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
                              ftor );

        // Need to fence before sending values off to be written out
        d_exec_space.fence();

        // copy rcols and vvals into std::vectors and write out
        out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    } else {
        d_exec_space.fence(); // still finish with a fence if no offd term present
    }
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::scale(
    AMP::Scalar alpha_in, MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple manually
    auto coeffs_d = std::get<3>( wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData ) );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();
    auto alpha        = static_cast<scalar_t>( alpha_in );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::scale (d)",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_d ),
        KOKKOS_LAMBDA( lidx_t n ) { coeffs_d( n ) *= alpha; } );

    if ( csrData->hasOffDiag() ) {
        auto coeffs_od =
            std::get<3>( wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData ) );

        const auto tnnz_od = csrData->numberOfNonZerosOffDiag();

        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::scale (od)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_od ),
            KOKKOS_LAMBDA( lidx_t n ) { coeffs_od( n ) *= alpha; } );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::matMultiply(
    MatrixData const &, MatrixData const &, MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsKokkos not implemented" );
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::axpy( AMP::Scalar alpha_in,
                                                                               const MatrixData &X,
                                                                               MatrixData &Y )
{
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( X ) );
    const auto csrDataY = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( Y ) );

    AMP_DEBUG_INSIST( csrDataX->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataY->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );
    AMP_DEBUG_INSIST( csrDataX->getMemoryLocation() == csrDataY->getMemoryLocation(),
                      "CSRMatrixOperationsKokkos::axpy X and Y must be in same memory space" );

    auto alpha = static_cast<scalar_t>( alpha_in );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto coeffsX_d =
            std::get<3>( wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrDataX ) );
        auto coeffsY_d =
            std::get<3>( wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrDataY ) );

        const auto tnnz_d = csrDataX->numberOfNonZerosDiag();
        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::axpy (d)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_d ),
            KOKKOS_LAMBDA( gidx_t n ) { coeffsY_d( n ) += alpha * coeffsX_d( n ); } );
    }

    if ( csrDataX->hasOffDiag() ) {
        const auto tnnz_od = csrDataX->numberOfNonZerosDiag();
        auto coeffsX_od =
            std::get<3>( wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrDataX ) );
        auto coeffsY_od =
            std::get<3>( wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrDataY ) );

        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::axpy (od)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_od ),
            KOKKOS_LAMBDA( gidx_t n ) { coeffsY_od( n ) += alpha * coeffsX_od( n ); } );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::setScalar(
    AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto alpha   = static_cast<scalar_t>( alpha_in );
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple manually
        auto coeffs_d =
            std::get<3>( wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData ) );
        Kokkos::deep_copy( d_exec_space, coeffs_d, alpha );
    }

    if ( csrData->hasOffDiag() ) {
        auto coeffs_od =
            std::get<3>( wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData ) );
        Kokkos::deep_copy( d_exec_space, coeffs_od, alpha );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::setDiagonal(
    std::shared_ptr<const Vector> in, MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow    = csrData->beginRow();

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple of views manually
    auto vtpl        = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
    auto nnz_d       = std::get<0>( vtpl );
    auto cols_d      = std::get<1>( vtpl );
    auto coeffs_d    = std::get<3>( vtpl );
    auto rowstarts_d = std::get<4>( vtpl );

    Kokkos::View<const scalar_t *, Kokkos::LayoutRight> vvals( in->getRawDataBlock<scalar_t>(),
                                                               nRows );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::setDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz_d( row );
            const auto rs = rowstarts_d( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols_d( rs + c ) == static_cast<gidx_t>( row + beginRow ) ) {
                    coeffs_d( rs + c ) = vvals( row );
                    break;
                }
            }
        } );

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::setIdentity(
    MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    zero( A );

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow    = csrData->beginRow();

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple of views manually
    auto vtpl        = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
    auto nnz_d       = std::get<0>( vtpl );
    auto cols_d      = std::get<1>( vtpl );
    auto coeffs_d    = std::get<3>( vtpl );
    auto rowstarts_d = std::get<4>( vtpl );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::setIdentity",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz_d( row );
            const auto rs = rowstarts_d( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols_d( rs + c ) == static_cast<gidx_t>( row + beginRow ) ) {
                    coeffs_d( rs + c ) = static_cast<scalar_t>( 1.0 );
                    break;
                }
            }
        } );

    d_exec_space.fence();
}


template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::extractDiagonal(
    MatrixData const &A, std::shared_ptr<Vector> buf )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow    = csrData->beginRow();

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple of views manually
    auto vtpl        = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
    auto nnz_d       = std::get<0>( vtpl );
    auto cols_d      = std::get<1>( vtpl );
    auto coeffs_d    = std::get<3>( vtpl );
    auto rowstarts_d = std::get<4>( vtpl );

    Kokkos::View<scalar_t *, Kokkos::LayoutRight> vvals( buf->getRawDataBlock<scalar_t>(), nRows );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::extractDiagonal",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
        KOKKOS_LAMBDA( lidx_t row ) {
            const auto nC = nnz_d( row );
            const auto rs = rowstarts_d( row );

            for ( lidx_t c = 0; c < nC; ++c ) {
                if ( cols_d( rs + c ) == static_cast<gidx_t>( row + beginRow ) ) {
                    vvals( row ) = coeffs_d( rs + c );
                    break;
                }
            }
        } );

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace, class ViewSpace>
AMP::Scalar CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace, ViewSpace>::LinfNorm(
    MatrixData const &A ) const
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    AMP_DEBUG_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                      "CSRMatrixOperationsKokkos is not implemented for device memory" );

    scalar_t rmax = 0.0;

    // Easier to write separate kernels for cases with(out) offd block
    if ( !csrData->hasOffDiag() ) {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl_d      = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_d       = std::get<0>( vtpl_d );
        auto coeffs_d    = std::get<3>( vtpl_d );
        auto rowstarts_d = std::get<4>( vtpl_d );

        Kokkos::parallel_reduce(
            "CSRMatrixOperationsKokkos::LinfNorm",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
            KOKKOS_LAMBDA( lidx_t row, scalar_t & lmax ) {
                scalar_t sum = 0;
                auto nC      = nnz_d( row );
                auto rs      = rowstarts_d( row );
                for ( lidx_t c = 0; c < nC; ++c ) {
                    sum += Kokkos::fabs( coeffs_d( rs + c ) );
                }
                lmax = lmax > sum ? lmax : sum;
            },
            Kokkos::Min<scalar_t>( rmax ) );
    } else {
        auto vtpl_d      = wrapCSRDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_d       = std::get<0>( vtpl_d );
        auto coeffs_d    = std::get<3>( vtpl_d );
        auto rowstarts_d = std::get<4>( vtpl_d );

        auto vtpl_od      = wrapCSROffDiagDataKokkos<Policy, Allocator, ViewSpace>( csrData );
        auto nnz_od       = std::get<0>( vtpl_od );
        auto coeffs_od    = std::get<3>( vtpl_od );
        auto rowstarts_od = std::get<4>( vtpl_od );

        Kokkos::parallel_reduce(
            "CSRMatrixOperationsKokkos::LinfNorm",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, nRows ),
            KOKKOS_LAMBDA( lidx_t row, scalar_t & lmax ) {
                scalar_t sum = 0;
                auto nC      = nnz_d( row );
                auto rs      = rowstarts_d( row );
                for ( lidx_t c = 0; c < nC; ++c ) {
                    sum += Kokkos::fabs( coeffs_d( rs + c ) );
                }
                nC = nnz_od( row );
                rs = rowstarts_od( row );
                for ( lidx_t c = 0; c < nC; ++c ) {
                    sum += Kokkos::fabs( coeffs_od( rs + c ) );
                }
                lmax = lmax > sum ? lmax : sum;
            },
            Kokkos::Min<scalar_t>( rmax ) );
    }

    // Reduce row sums to get global Linf norm
    AMP_MPI comm = csrData->getComm();
    return comm.maxReduce<scalar_t>( rmax );
}

} // namespace AMP::LinearAlgebra

#endif // close check for Kokkos being defined
