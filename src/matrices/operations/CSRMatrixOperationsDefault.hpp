#include "AMP/matrices/data/CSRMatrixData.h"
#include "AMP/matrices/operations/CSRMatrixOperationsDefault.h"
#include "AMP/utils/Utilities.h"

#include <algorithm>

namespace AMP::LinearAlgebra {

template<typename Policy>
static CSRMatrixData<Policy> const *getCSRMatrixData( MatrixData const &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> const *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}

template<typename Policy>
static CSRMatrixData<Policy> *getCSRMatrixData( MatrixData &A )
{
    auto ptr = dynamic_cast<CSRMatrixData<Policy> *>( &A );
    AMP_INSIST( ptr, "dynamic cast from const MatrixData to const CSRMatrixData failed" );
    return ptr;
}


template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::mult( std::shared_ptr<const Vector> in,
                                               MatrixData const &A,
                                               std::shared_ptr<Vector> out )
{
    AMP_ASSERT( in && out );

    using gidx_t   = typename Policy::gidx_t;
    using lidx_t   = typename Policy::lidx_t;
    using scalar_t = typename Policy::scalar_t;

    auto csrData = getCSRMatrixData<Policy>( const_cast<MatrixData &>( A ) );

    auto [nnz, cols, coeffs] = csrData->getCSRData();

    auto memType = AMP::Utilities::getMemoryType( cols );
    AMP_INSIST( memType == AMP::Utilities::MemoryType::host ||
                    memType == AMP::Utilities::MemoryType::unregistered,
                "CSRMatrixOperationsDefault is implemented only for host memory" );

    const auto nRows = csrData->numLocalRows();
    auto maxColLen   = *std::max_element( nnz, nnz + nRows );

    std::vector<size_t> rcols( maxColLen );
    std::vector<scalar_t> vvals( maxColLen );

    auto beginRow = csrData->beginRow();

    lidx_t offset = 0;
    for ( lidx_t row = 0; row < nRows; ++row ) {

        const auto nCols = nnz[row];

        const auto cloc = &cols[offset];
        const auto vloc = &coeffs[offset];

        std::transform(
            cloc, cloc + nCols, rcols.begin(), []( gidx_t col ) -> size_t { return col; } );

        in->getValuesByGlobalID( nCols, rcols.data(), vvals.data() );

        scalar_t val =
            std::inner_product( vloc, vloc + nCols, vvals.data(), static_cast<scalar_t>( 0.0 ) );

        out->setValueByGlobalID( static_cast<size_t>( beginRow + row ), val );

        offset += nCols;
    }

    out->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::multTranspose( std::shared_ptr<const Vector> in,
                                                        MatrixData const &A,
                                                        std::shared_ptr<Vector> out )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::scale( AMP::Scalar alpha_in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::matMultiply( MatrixData const &Am,
                                                      MatrixData const &Bm,
                                                      MatrixData &Cm )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::axpy( AMP::Scalar alpha_in,
                                               const MatrixData &X,
                                               MatrixData &Y )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setScalar( AMP::Scalar alpha_in, MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::zero( MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setDiagonal( std::shared_ptr<const Vector> in,
                                                      MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
void CSRMatrixOperationsDefault<Policy>::setIdentity( MatrixData &A )
{
    AMP_ERROR( "Not implemented" );
}

template<typename Policy>
AMP::Scalar CSRMatrixOperationsDefault<Policy>::L1Norm( MatrixData const &A ) const
{
    AMP_ERROR( "Not implemented" );
}

} // namespace AMP::LinearAlgebra
