#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined(AMP_USE_KOKKOS) || defined(AMP_USE_TRILINOS_KOKKOS)

#include "Kokkos_core.h"

namespace AMP::LinearAlgebra {

template<typename Policy>
auto wrapCSRDiagDataKokkos( std::shared_ptr<CSRMatrixData<Policy>> csrData )
{
  using lidx_t   = typename Policy::lidx_t;
  using gidx_t   = typename Policy::gidx_t;
  using scalar_t = typename Policy::scalar_t;

  const lidx_t nrows = static_cast<lidx_t>( csrData->numLocalRows() );
  const lidx_t nnz_tot = csrData->numberOfNonZerosDiag();
  auto [nnz, cols, cols_loc, coeffs] = csrData->getCSRDiagData();

  return std::make_tuple( Kokkos::View<lidx_t*,Kokkos::LayoutRight>( nnz, nrows ),
			  Kokkos::View<gidx_t*,Kokkos::LayoutRight>( cols, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( cols_loc, nnz_tot ),
			  Kokkos::View<scalar_t*,Kokkos::LayoutRight>( coeffs, nnz_tot ) );
}
  
template<typename Policy>
auto wrapCSROffDiagDataKokkos( std::shared_ptr<CSRMatrixData<Policy>> csrData )
{
  using lidx_t   = typename Policy::lidx_t;
  using gidx_t   = typename Policy::gidx_t;
  using scalar_t = typename Policy::scalar_t;

  const lidx_t nrows = static_cast<lidx_t>( csrData->numLocalRows() );
  const lidx_t nnz_tot = csrData->numberOfNonZerosOffDiag();
  auto [nnz, cols, cols_loc, coeffs] = csrData->getCSROffDiagData();

  return std::make_tuple( Kokkos::View<lidx_t*,Kokkos::LayoutRight>( nnz, nrows ),
			  Kokkos::View<gidx_t*,Kokkos::LayoutRight>( cols, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( cols_loc, nnz_tot ),
			  Kokkos::View<scalar_t*,Kokkos::LayoutRight>( coeffs, nnz_tot ) );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::mult( std::shared_ptr<const Vector> in,
					      MatrixData const &A,
					      std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsKokkos::mult" );
    AMP_DEBUG_ASSERT( in && out );
    AMP_DEBUG_ASSERT( in->getUpdateStatus() == AMP::LinearAlgebra::UpdateState::UNCHANGED );

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = wrapCSRDiagDataKokkos( csrData );

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );
    const auto &ghosts          = inData->getGhosts();
    auto outData                = out->getVectorData();
    scalar_t *outDataBlock      = outData->getRawDataBlock<scalar_t>( 0 );

    AMP_INSIST( AMP::Utilities::getMemoryType( cols_loc_d ) != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    AMP_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with one data block" );

    AMP_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with accessible ghosts" );

    AMP_ASSERT( inDataBlock && outDataBlock );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    {
        PROFILE( "CSRMatrixOperationsKokkos::mult (local)" );
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {
            const auto nCols = nnz_d[row];
            const auto cloc  = &cols_loc_d[offset];
            const auto vloc  = &coeffs_d[offset];

            outDataBlock[row] = 0.0;
            for ( lidx_t c = 0; c < nCols; ++c ) {
                outDataBlock[row] += vloc[c] * inDataBlock[cloc[c]];
            }

            offset += nCols;
        }
    }

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsKokkos::mult (ghost)" );
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        lidx_t offset                                  = 0;

        for ( lidx_t row = 0; row < nRows; ++row ) {
            const auto nCols = nnz_od[row];
            const auto cloc  = &cols_loc_od[offset];
            const auto vloc  = &coeffs_od[offset];

            for ( lidx_t c = 0; c < nCols; ++c ) {
                outDataBlock[row] += vloc[c] * ghosts[cloc[c]];
            }

            offset += nCols;
        }
    }
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::multTranspose( std::shared_ptr<const Vector> in,
						       MatrixData const &A,
						       std::shared_ptr<Vector> out )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::matMultiply( MatrixData const &,
						     MatrixData const &,
						     MatrixData & )
{
    AMP_WARNING( "SpGEMM for CSRMatrixOperationsKokkos not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::axpy( AMP::Scalar alpha_in,
                                               const MatrixData &X,
                                               MatrixData &Y )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d] = csrData->getCSRDiagData();

    auto memType = AMP::Utilities::getMemoryType( cols_loc_d );
    AMP_INSIST( memType != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();

    auto alpha = static_cast<scalar_t>( alpha_in );
    std::fill(
        const_cast<scalar_t *>( coeffs_d ), const_cast<scalar_t *>( coeffs_d ) + tnnz_d, alpha );

    if ( csrData->hasOffDiag() ) {
        const auto tnnz_od                             = csrData->numberOfNonZerosOffDiag();
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od] = csrData->getCSROffDiagData();
        std::fill( const_cast<scalar_t *>( coeffs_od ),
                   const_cast<scalar_t *>( coeffs_od ) + tnnz_od,
                   alpha );
    }
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::zero( MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;
    setScalar( static_cast<scalar_t>( 0.0 ), A );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::setDiagonal( std::shared_ptr<const Vector> in,
						     MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::setIdentity( MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
AMP::Scalar CSRMatrixOperationsKokkos<Policy>::L1Norm( MatrixData const &A ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra

#endif // close check for Kokkos being defined
