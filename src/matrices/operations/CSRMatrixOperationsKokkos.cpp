#include "AMP/AMP_TPLs.h"
#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.h"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkosFunctors.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined( AMP_USE_KOKKOS ) || defined( AMP_USE_TRILINOS_KOKKOS )

    #include "Kokkos_Core.hpp"

namespace AMP::LinearAlgebra {

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::mult(
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
    // const auto nnz_per_row_d = csrData->numberOfNonZerosDiag() / nRows;
    // const auto nnz_per_row_od = csrData->numberOfNonZerosOffDiag() / nRows;

    auto inData        = in->getVectorData();
    const auto &ghosts = inData->getGhosts();
    const auto nGhosts = static_cast<lidx_t>( ghosts.size() );
    auto outData       = out->getVectorData();

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );

    AMP_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with one data block" );

    AMP_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with accessible ghosts" );


    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 Kokkos::SharedSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        inDataBlock( inData->getRawDataBlock<scalar_t>( 0 ), nCols );
    Kokkos::View<scalar_t *, Kokkos::LayoutRight, Kokkos::SharedSpace> outDataBlock(
        outData->getRawDataBlock<scalar_t>( 0 ), nRows );
    Kokkos::View<const scalar_t *,
                 Kokkos::LayoutRight,
                 Kokkos::SharedSpace,
                 Kokkos::MemoryTraits<Kokkos::RandomAccess>>
        ghostDataBlock( ghosts.data(), nGhosts );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl        = wrapCSRDiagDataKokkos( csrData );
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

        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> team_policy(
            d_exec_space, num_teams, Kokkos::AUTO, vector_length );

        Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (local)", team_policy, ftor );
    }

    if ( csrData->hasOffDiag() ) {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl         = wrapCSROffDiagDataKokkos( csrData );
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

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::multTranspose(
    std::shared_ptr<const Vector> in, MatrixData const &A, std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsKokkos::multTranspose" );

    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    const auto nCols = static_cast<lidx_t>( csrData->numLocalColumns() );
    auto inData      = in->getVectorData();
    auto outData     = out->getVectorData();

    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t *, Kokkos::LayoutRight, Kokkos::SharedSpace> inDataBlock(
        inData->getRawDataBlock<scalar_t>( 0 ), nRows );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl        = wrapCSRDiagDataKokkos( csrData );
        auto nnz_d       = std::get<0>( vtpl );
        auto cols_loc_d  = std::get<2>( vtpl );
        auto coeffs_d    = std::get<3>( vtpl );
        auto rowstarts_d = std::get<4>( vtpl );

        // Make temporary views for output columns and values
        auto outDataBlock = Kokkos::View<scalar_t *,
                                         Kokkos::LayoutRight,
                                         Kokkos::SharedSpace,
                                         Kokkos::MemoryTraits<Kokkos::Atomic>>(
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

        Kokkos::TeamPolicy<ExecSpace, Kokkos::Schedule<Kokkos::Dynamic>> team_policy(
            d_exec_space, num_teams, Kokkos::AUTO, vector_length );

        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::multTranspose (local)", team_policy, ftor );
    }

    if ( csrData->hasOffDiag() ) {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto vtpl         = wrapCSROffDiagDataKokkos( csrData );
        auto nnz_od       = std::get<0>( vtpl );
        auto cols_loc_od  = std::get<2>( vtpl );
        auto coeffs_od    = std::get<3>( vtpl );
        auto rowstarts_od = std::get<4>( vtpl );

        // get diag map but leave in std::vector
        // it is not needed inside compute kernel this time
        std::vector<size_t> rcols;
        csrData->getOffDiagColumnMap( rcols );
        const auto num_unq = rcols.size();

        // Make temporary view for output values
        auto vvals = Kokkos::View<scalar_t *,
                                  Kokkos::LayoutRight,
                                  Kokkos::SharedSpace,
                                  Kokkos::MemoryTraits<Kokkos::Atomic>>( "vvals", num_unq );
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

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::scale( AMP::Scalar alpha_in,
                                                                     MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple manually
    auto coeffs_d = std::get<3>( wrapCSRDiagDataKokkos( csrData ) );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();
    auto alpha        = static_cast<scalar_t>( alpha_in );

    Kokkos::parallel_for(
        "CSRMatrixOperationsKokkos::scale (d)",
        Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_d ),
        KOKKOS_LAMBDA( lidx_t n ) { coeffs_d( n ) *= alpha; } );

    if ( csrData->hasOffDiag() ) {
        auto coeffs_od = std::get<3>( wrapCSROffDiagDataKokkos( csrData ) );

        const auto tnnz_od = csrData->numberOfNonZerosOffDiag();

        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::scale (od)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_od ),
            KOKKOS_LAMBDA( lidx_t n ) { coeffs_od( n ) *= alpha; } );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::matMultiply( MatrixData const &,
                                                                           MatrixData const &,
                                                                           MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsKokkos not implemented" );
}

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::axpy( AMP::Scalar alpha_in,
                                                                    const MatrixData &X,
                                                                    MatrixData &Y )
{
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( X ) );
    const auto csrDataY = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( Y ) );

    AMP_INSIST( csrDataX->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );
    AMP_INSIST( csrDataY->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );
    AMP_INSIST( csrDataX->getMemoryLocation() == csrDataY->getMemoryLocation(),
                "CSRMatrixOperationsKokkos::axpy X and Y must be in same memory space" );

    auto alpha = static_cast<scalar_t>( alpha_in );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple of views manually
        auto coeffsX_d = std::get<3>( wrapCSRDiagDataKokkos( csrDataX ) );
        auto coeffsY_d = std::get<3>( wrapCSRDiagDataKokkos( csrDataY ) );

        const auto tnnz_d = csrDataX->numberOfNonZerosDiag();
        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::axpy (d)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_d ),
            KOKKOS_LAMBDA( gidx_t n ) { coeffsY_d( n ) += alpha * coeffsX_d( n ); } );
    }

    if ( csrDataX->hasOffDiag() ) {
        const auto tnnz_od = csrDataX->numberOfNonZerosDiag();
        auto coeffsX_od    = std::get<3>( wrapCSROffDiagDataKokkos( csrDataX ) );
        auto coeffsY_od    = std::get<3>( wrapCSROffDiagDataKokkos( csrDataY ) );

        Kokkos::parallel_for(
            "CSRMatrixOperationsKokkos::axpy (od)",
            Kokkos::RangePolicy<ExecSpace>( d_exec_space, 0, tnnz_od ),
            KOKKOS_LAMBDA( gidx_t n ) { coeffsY_od( n ) += alpha * coeffsX_od( n ); } );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::setScalar( AMP::Scalar alpha_in,
                                                                         MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto alpha   = static_cast<scalar_t>( alpha_in );
    auto csrData = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );

    {
        // lambda capture of structured bindings throws warning on c++17
        // unpack the tuple manually
        auto coeffs_d = std::get<3>( wrapCSRDiagDataKokkos( csrData ) );
        Kokkos::deep_copy( d_exec_space, coeffs_d, alpha );
    }

    if ( csrData->hasOffDiag() ) {
        auto coeffs_od = std::get<3>( wrapCSROffDiagDataKokkos( csrData ) );
        Kokkos::deep_copy( d_exec_space, coeffs_od, alpha );
    }

    d_exec_space.fence();
}

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::setDiagonal(
    std::shared_ptr<const Vector> in, MatrixData &A )
{
    using lidx_t   = typename Policy::lidx_t;
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData     = getCSRMatrixData<Policy, Allocator>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow    = csrData->beginRow();

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is not implemented for device memory" );

    // lambda capture of structured bindings throws warning on c++17
    // unpack the tuple of views manually
    auto vtpl        = wrapCSRDiagDataKokkos( csrData );
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

template<typename Policy, typename Allocator, class ExecSpace>
void CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::setIdentity( MatrixData &A )
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
    auto vtpl        = wrapCSRDiagDataKokkos( csrData );
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

template<typename Policy, typename Allocator, class ExecSpace>
AMP::Scalar
CSRMatrixOperationsKokkos<Policy, Allocator, ExecSpace>::L1Norm( MatrixData const & ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra

#endif // close check for Kokkos being defined
