#include "ProfilerApp.h"

namespace AMP::LinearAlgebra {

template<class BIGINT_TYPE, class INT_TYPE>
void GetRowHelper::NNZ( BIGINT_TYPE row, INT_TYPE &num_local, INT_TYPE &num_remote )
{
    PROFILE( "GetRowHelper::NNZ" );

    // attempt to get row dofs and resize backing vector if needed
    const auto id = d_leftDOF->getElementID( row );
    const auto N  = d_rightDOF->getRowDOFs( id, d_rowDOFs.data(), d_rowDOFs.size(), false );
    if ( N > d_rowDOFs.size() ) {
        d_rowDOFs.resize( N );
        d_rightDOF->getRowDOFs( id, d_rowDOFs.data(), d_rowDOFs.size(), false );
    }

    // test for local vs remote and assign accordingly
    num_local  = 0;
    num_remote = 0;
    for ( size_t k = 0; k < N; ++k ) {
        if ( d_rowDOFs[k] >= d_beginCol && d_rowDOFs[k] < d_endCol ) {
            ++num_local;
        } else {
            ++num_remote;
        }
    }
}

template<class BIGINT_TYPE>
void GetRowHelper::getRow( BIGINT_TYPE row, BIGINT_TYPE *cols_local, BIGINT_TYPE *cols_remote )
{
    PROFILE( "GetRowHelper::getRow" );

    // attempt to get row dofs and resize backing vector if needed
    // resize shouldn't be needed since a prior call to NNZ should have already
    // done this
    const auto id = d_leftDOF->getElementID( row );
    const auto N  = d_rightDOF->getRowDOFs( id, d_rowDOFs.data(), d_rowDOFs.size(), false );
    if ( N > d_rowDOFs.size() ) {
        AMP_WARNING(
            "GetRowHelper: resize of backing vector needed in getRow. This should not happen." );
        d_rowDOFs.resize( N );
        d_rightDOF->getRowDOFs( id, d_rowDOFs.data(), d_rowDOFs.size(), false );
    }

    // test for local vs remote and assign accordingly
    size_t nl = 0, nr = 0;
    for ( size_t k = 0; k < N; ++k ) {
        if ( d_rowDOFs[k] >= d_beginCol && d_rowDOFs[k] < d_endCol ) {
            cols_local[nl++] = static_cast<BIGINT_TYPE>( d_rowDOFs[k] );
        } else {
            cols_remote[nr++] = static_cast<BIGINT_TYPE>( d_rowDOFs[k] );
        }
    }
}

} // namespace AMP::LinearAlgebra
