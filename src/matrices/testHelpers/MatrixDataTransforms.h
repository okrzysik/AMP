#ifndef included_MatrixDataTransforms_h
#define included_MatrixDataTransforms_h

#include "AMP/matrices/Matrix.h"

#include <type_traits>

namespace AMP::LinearAlgebra {

template<typename Policy>
void transformDofToCSR( std::shared_ptr<Matrix> matrix,
                        typename Policy::gidx_t &firstRow,
                        typename Policy::gidx_t &endRow,
                        std::vector<typename Policy::lidx_t> &nnz,
                        std::vector<typename Policy::gidx_t> &cols,
                        std::vector<typename Policy::scalar_t> &coeffs )
{
    AMP_ASSERT( matrix );

    firstRow = static_cast<typename Policy::gidx_t>( matrix->beginRow() );
    endRow   = static_cast<typename Policy::gidx_t>( matrix->endRow() );

    for ( auto row = firstRow; row < endRow; ++row ) {

        std::vector<size_t> rcols;
        std::vector<typename Policy::scalar_t> rvals;

        matrix->getRowByGlobalID( row, rcols, rvals );
        nnz.push_back( static_cast<typename Policy::lidx_t>( rcols.size() ) );

        if ( std::is_same_v<size_t, typename Policy::gidx_t> ) {
            cols.insert( cols.end(), rcols.begin(), rcols.end() );
        } else {
            std::transform( rcols.cbegin(),
                            rcols.cend(),
                            std::back_inserter( cols ),
                            []( size_t col ) -> typename Policy::gidx_t { return col; } );
        }

        coeffs.insert( coeffs.end(), rvals.cbegin(), rvals.cend() );
    }
}
} // namespace AMP::LinearAlgebra

#endif
