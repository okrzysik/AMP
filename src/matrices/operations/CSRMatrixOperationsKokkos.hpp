#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsKokkos.h"
#include "AMP/utils/Utilities.h"
#include "AMP/AMP_TPLs.h"

#include <algorithm>

#include "ProfilerApp.h"

#if defined(AMP_USE_KOKKOS) || defined(AMP_USE_TRILINOS_KOKKOS)

#include "Kokkos_Core.hpp"

namespace AMP::LinearAlgebra {

template<typename Policy>
auto wrapCSRDiagDataKokkos( CSRMatrixData<Policy> *csrData )
{
  using lidx_t   = typename Policy::lidx_t;
  using gidx_t   = typename Policy::gidx_t;
  using scalar_t = typename Policy::scalar_t;

  const lidx_t nrows = static_cast<lidx_t>( csrData->numLocalRows() );
  const lidx_t nnz_tot = csrData->numberOfNonZerosDiag();
  auto [nnz, cols, cols_loc, coeffs] = csrData->getCSRDiagData();
  auto *rowstarts = csrData->getDiagRowStarts();

  return std::make_tuple( Kokkos::View<lidx_t*,Kokkos::LayoutRight>( nnz, nrows ),
			  Kokkos::View<gidx_t*,Kokkos::LayoutRight>( cols, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( cols_loc, nnz_tot ),
			  Kokkos::View<scalar_t*,Kokkos::LayoutRight>( coeffs, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( rowstarts, nrows ) );
}
  
template<typename Policy>
auto wrapCSROffDiagDataKokkos( CSRMatrixData<Policy> *csrData )
{
  using lidx_t   = typename Policy::lidx_t;
  using gidx_t   = typename Policy::gidx_t;
  using scalar_t = typename Policy::scalar_t;

  const lidx_t nrows = static_cast<lidx_t>( csrData->numLocalRows() );
  const lidx_t nnz_tot = csrData->numberOfNonZerosOffDiag();
  auto [nnz, cols, cols_loc, coeffs] = csrData->getCSROffDiagData();
  auto *rowstarts = csrData->getOffDiagRowStarts();

  return std::make_tuple( Kokkos::View<lidx_t*,Kokkos::LayoutRight>( nnz, nrows ),
			  Kokkos::View<gidx_t*,Kokkos::LayoutRight>( cols, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( cols_loc, nnz_tot ),
			  Kokkos::View<scalar_t*,Kokkos::LayoutRight>( coeffs, nnz_tot ),
			  Kokkos::View<lidx_t*,Kokkos::LayoutRight>( rowstarts, nrows ) );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::mult( std::shared_ptr<const Vector> in,
					      MatrixData const &A,
					      std::shared_ptr<Vector> out )
{
    PROFILE( "CSRMatrixOperationsKokkos::mult" );
    AMP_DEBUG_ASSERT( in && out );

    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    const auto nCols = static_cast<lidx_t>( csrData->numLocalColumns() );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d, rowstarts_d] = wrapCSRDiagDataKokkos( csrData );

    auto inData                 = in->getVectorData();
    const auto &ghosts          = inData->getGhosts();
    const auto nGhosts          = static_cast<lidx_t>( ghosts.size() );
    auto outData                = out->getVectorData();

    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    AMP_INSIST(
        1 == inData->numberOfDataBlocks(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with one data block" );

    AMP_INSIST(
        ghosts.size() == inData->getGhostSize(),
        "CSRMatrixOperationsKokkos::mult only implemented for vectors with accessible ghosts" );
    
    // Wrap in/out data into Kokkos Views
    Kokkos::View<const scalar_t*,Kokkos::LayoutRight> inDataBlock( inData->getRawDataBlock<scalar_t>( 0 ), nCols );
    Kokkos::View<scalar_t*,Kokkos::LayoutRight> outDataBlock( outData->getRawDataBlock<scalar_t>( 0 ), nRows );
    Kokkos::View<const scalar_t*,Kokkos::LayoutRight> ghostDataBlock( ghosts.data(), nGhosts );

    {
        PROFILE( "CSRMatrixOperationsKokkos::mult (local)" );

	Kokkos::parallel_for( Kokkos::TeamPolicy<Kokkos::IndexType<lidx_t>>(nRows, 1, 4),
			      KOKKOS_LAMBDA (const Kokkos::TeamPolicy<>::member_type &tm) {
				lidx_t row = tm.league_rank();
				const auto nC = nnz_d( row );
				const auto rs = rowstarts_d( row );
				scalar_t sum = 0.0;
				Kokkos::parallel_reduce( Kokkos::TeamVectorRange(tm,nC),
							 [=] (lidx_t &c, scalar_t &lsum) {
							   const auto cl = cols_loc_d( rs + c );
							   lsum += coeffs_d( rs + c ) * inDataBlock( cl );
							 }, sum);
				outDataBlock(row) = sum;
			      } );
    }

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsKokkos::mult (ghost)" );
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od, rowstarts_od] = wrapCSROffDiagDataKokkos( csrData );

	Kokkos::parallel_for( "CSRMatrixOperationsKokkos::mult (ghost)",
			      nRows, KOKKOS_LAMBDA ( const lidx_t row ) {
				const auto nC = nnz_od( row );
				const auto rs = rowstarts_od( row );
				
				for ( lidx_t c = 0; c < nC; ++c ) {
				  const auto cl = cols_loc_od( rs + c );
				  outDataBlock( row ) += coeffs_od( rs + c ) * ghostDataBlock( cl );
				}
			      } );
    }
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::multTranspose( std::shared_ptr<const Vector> in,
						       MatrixData const &A,
						       std::shared_ptr<Vector> out)
{
    PROFILE( "CSRMatrixOperationsKokkos::multTranspose" );
    
    // this is not meant to be an optimized version. It is provided for completeness
    AMP_DEBUG_ASSERT( in && out );

    out->zero();

    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    
    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    auto inData                 = in->getVectorData();
    const scalar_t *inDataBlock = inData->getRawDataBlock<scalar_t>( 0 );

    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );

    {
        PROFILE( "CSRMatrixOperationsKokkos::multTranspose (d)" );
	auto [nnz, cols, cols_loc, coeffs] = csrData->getCSRDiagData();
	auto [num_unq, cols_unq] = csrData->getDiagColumnMap();
	
        std::vector<scalar_t> vvals( num_unq, 0.0 );
	
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto ncols = nnz[row];
            const auto cloc = &cols[offset];
            const auto vloc = &coeffs[offset];
            const auto val = inDataBlock[ row ];

            for ( lidx_t j = 0; j < ncols; ++j ) {
	        auto cit = std::lower_bound( cols_unq, cols_unq + num_unq, cloc[j] );
		auto idx = std::distance( cols_unq, cit );
                vvals[idx] += vloc[j] * val;
            }

            offset += ncols;
        }

	// convert colmap to size_t and write out data
	std::vector<size_t> rcols( num_unq );
	std::transform( cols_unq, cols_unq + num_unq, rcols.begin(),
			[]( gidx_t col ) -> size_t { return col; } );
	out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
    out->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );

    if ( csrData->hasOffDiag() ) {
        PROFILE( "CSRMatrixOperationsKokkos::multTranspose (od)" );
	auto [nnz, cols, cols_loc, coeffs] = csrData->getCSROffDiagData();
	auto [num_unq, cols_unq] = csrData->getOffDiagColumnMap();
	
        std::vector<scalar_t> vvals( num_unq, 0.0 );
	
        lidx_t offset = 0;
        for ( lidx_t row = 0; row < nRows; ++row ) {

            const auto ncols = nnz[row];
	    if ( ncols == 0 ) { continue; }
	    
            const auto cloc = &cols[offset];
            const auto vloc = &coeffs[offset];
            const auto val = inDataBlock[ row ];

            for ( lidx_t j = 0; j < ncols; ++j ) {
	        auto cit = std::lower_bound( cols_unq, cols_unq + num_unq, cloc[j] );
		auto idx = std::distance( cols_unq, cit );
                vvals[idx] += vloc[j] * val;
            }

            offset += ncols;
        }

	// convert colmap to size_t and write out data
	std::vector<size_t> rcols( num_unq );
	std::transform( cols_unq, cols_unq + num_unq, rcols.begin(),
			[]( gidx_t col ) -> size_t { return col; } );
	out->addValuesByGlobalID( num_unq, rcols.data(), vvals.data() );
    }
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    using lidx_t = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;
    
    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    
    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d, rowstarts_d] = wrapCSRDiagDataKokkos( csrData );

    const auto tnnz_d = csrData->numberOfNonZerosDiag();
    auto alpha = static_cast<scalar_t>( alpha_in );
    
    Kokkos::parallel_for( "CSRMatrixOperationsKokkos::scale (d)",
			  tnnz_d, KOKKOS_LAMBDA ( const lidx_t n ) {
			    coeffs_d( n ) *= alpha;
			  } );

    if ( csrData->hasOffDiag() ) {
    auto [nnz_od, cols_od, cols_loc_od, coeffs_od, rowstarts_od] = wrapCSROffDiagDataKokkos( csrData );

    const auto tnnz_od = csrData->numberOfNonZerosOffDiag();
    
    Kokkos::parallel_for( "CSRMatrixOperationsKokkos::scale (od)",
			  tnnz_od, KOKKOS_LAMBDA ( const lidx_t n ) {
			    coeffs_od( n ) *= alpha;
			  } );
    }
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
    using gidx_t   = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;

    const auto csrDataX = getCSRMatrixData<Policy>( const_cast<MatrixData &>( X ) );
    const auto csrDataY = getCSRMatrixData<Policy>( const_cast<MatrixData &>( Y ) );
    
    AMP_INSIST( csrDataX->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );
    AMP_INSIST( csrDataY->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );
    AMP_INSIST( csrDataX->getMemoryLocation() == csrDataY->getMemoryLocation(),
                "CSRMatrixOperationsKokkos::axpy X and Y must be in same memory space" );

    auto [nnzX_d, colsX_d, cols_locX_d, coeffsX_d, rowstartsX_d] = wrapCSRDiagDataKokkos( csrDataX );
    auto [nnzY_d, colsY_d, cols_locY_d, coeffsY_d, rowstartsY_d] = wrapCSRDiagDataKokkos( csrDataY );

    auto alpha = static_cast<scalar_t>( alpha_in );

    const auto tnnz_d = csrDataX->numberOfNonZerosDiag();
    Kokkos::parallel_for( "CSRMatrixOperationsKokkos::axpy (d)",
			  tnnz_d, KOKKOS_LAMBDA ( const gidx_t n ) {
			    coeffsY_d( n ) += alpha * coeffsX_d( n );
			  } );

    if ( csrDataX->hasOffDiag() ) {
      const auto tnnz_od = csrDataX->numberOfNonZerosDiag();
      auto [nnzX_od, colsX_od, cols_locX_od, coeffsX_od, rowstartsX_od] =
	wrapCSROffDiagDataKokkos( csrDataX );
      auto [nnzY_od, colsY_od, cols_locY_od, coeffsY_od, rowstartsY_od] =
	wrapCSROffDiagDataKokkos( csrDataY );
      
      Kokkos::parallel_for( "CSRMatrixOperationsKokkos::axpy (od)",
			    tnnz_od, KOKKOS_LAMBDA ( const gidx_t n ) {
			      coeffsY_od( n ) += alpha * coeffsX_od( n );
			    } );
    }
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    using scalar_t = typename Policy::scalar_t;

    auto alpha = static_cast<scalar_t>( alpha_in );
    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    
    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d, rowstarts_d] = wrapCSRDiagDataKokkos( csrData );

    Kokkos::deep_copy( coeffs_d, alpha );

    if ( csrData->hasOffDiag() ) {
        auto [nnz_od, cols_od, cols_loc_od, coeffs_od, rowstarts_od] = wrapCSROffDiagDataKokkos( csrData );
	Kokkos::deep_copy( coeffs_od, alpha );
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
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;
    
    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow = csrData->beginRow();
    
    AMP_INSIST( csrData->getMemoryLocation() != AMP::Utilities::MemoryType::device,
                "CSRMatrixOperationsKokkos is implemented only for host memory" );

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d, rowstarts_d] = wrapCSRDiagDataKokkos( csrData );

    Kokkos::View<const scalar_t*,Kokkos::LayoutRight> vvals( in->getRawDataBlock<scalar_t>(), nRows );

    Kokkos::parallel_for( "CSRMatrixOperationsKokkos::setDiagonal",
			  nRows, KOKKOS_LAMBDA ( const lidx_t row ) {
			    const auto nC = nnz_d( row );
			    const auto rs = rowstarts_d( row );
			    
			    for ( lidx_t c = 0; c < nC; ++c ) {
			      if ( cols_d( rs + c ) == static_cast<gidx_t>( row + beginRow ) ) {
				coeffs_d( rs + c ) = vvals( row );
				break;
			      }
			    }
			  } );
}

template<typename Policy>
void CSRMatrixOperationsKokkos<Policy>::setIdentity( MatrixData &A )
{
    using lidx_t = typename Policy::lidx_t;
    using gidx_t = typename Policy::gidx_t;
    using scalar_t = typename Policy::scalar_t;
    
    zero( A );
    
    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );
    const auto nRows = static_cast<lidx_t>( csrData->numLocalRows() );
    auto beginRow = csrData->beginRow();

    auto [nnz_d, cols_d, cols_loc_d, coeffs_d, rowstarts_d] = wrapCSRDiagDataKokkos( csrData );

    Kokkos::parallel_for( "CSRMatrixOperationsKokkos::setIdentity",
			  nRows, KOKKOS_LAMBDA ( const lidx_t row ) {
			    const auto nC = nnz_d( row );
			    const auto rs = rowstarts_d( row );
			    
			    for ( lidx_t c = 0; c < nC; ++c ) {
			      if ( cols_d( rs + c ) == static_cast<gidx_t>( row + beginRow ) ) {
				coeffs_d( rs + c ) = static_cast<scalar_t>( 1.0 );
				break;
			      }
			    }
			  } );
}

template<typename Policy>
AMP::Scalar CSRMatrixOperationsKokkos<Policy>::L1Norm( MatrixData const & ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra

#endif // close check for Kokkos being defined
