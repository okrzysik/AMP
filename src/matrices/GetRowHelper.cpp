#include "AMP/matrices/GetRowHelper.h"

#include <cstdlib>


namespace AMP::LinearAlgebra {


/****************************************************************
 * Constructor                                                   *
 ****************************************************************/
GetRowHelper::GetRowHelper( std::shared_ptr<const AMP::Discretization::DOFManager> leftDOF,
                            std::shared_ptr<const AMP::Discretization::DOFManager> rightDOF )
    : d_leftDOF( leftDOF ), d_rightDOF( rightDOF )
{
    size_t beginRow = d_leftDOF->beginDOF();
    size_t endRow   = d_leftDOF->endDOF();
    size_t N_rows   = endRow - beginRow;
    d_NNZ           = new std::array<size_t, 2>[N_rows];
    size_t beginCol = d_rightDOF->beginDOF();
    size_t endCol   = d_rightDOF->endDOF();
    std::vector<size_t> dofs( 1024 );
    for ( size_t row = beginRow, i = 0; i < N_rows; i++, row++ ) {
        d_NNZ[i] = { 0, 0 };
        auto id  = d_leftDOF->getElementID( row );
        size_t N = d_rightDOF->getRowDOFs( id, dofs.data(), dofs.size() );
        if ( N > dofs.size() ) {
            dofs.resize( N );
            N = d_rightDOF->getRowDOFs( id, dofs.data(), dofs.size() );
        }
        reserve( N );
        for ( size_t k = 0; k < N; k++ ) {
            if ( dofs[k] >= beginCol && dofs[k] < endCol ) {
                d_NNZ[i][0]++;
                d_local[d_size[0]++] = dofs[k];
            } else {
                d_NNZ[i][1]++;
                d_remote[d_size[1]++] = dofs[k];
            }
        }
    }
    d_localOffset     = new size_t[N_rows];
    d_remoteOffset    = new size_t[N_rows];
    d_localOffset[0]  = 0;
    d_remoteOffset[0] = 0;
    for ( size_t i = 1; i < N_rows; i++ ) {
        d_localOffset[i]  = d_localOffset[i - 1] + d_NNZ[i - 1][0];
        d_remoteOffset[i] = d_remoteOffset[i - 1] + d_NNZ[i - 1][1];
    }
}
GetRowHelper::~GetRowHelper()
{
    delete[] d_NNZ;
    delete[] d_localOffset;
    delete[] d_remoteOffset;
    free( d_local );
    free( d_remote );
}


/****************************************************************
 * Reserve additional capacity for dofs                          *
 ****************************************************************/
void GetRowHelper::reserve( size_t N )
{
    if ( d_size[0] + N > d_capacity[0] ) {
        while ( d_size[0] + N > d_capacity[0] )
            d_capacity[0] = std::max<size_t>( 1024, 2 * d_capacity[0] );
        d_local = (size_t *) realloc( d_local, d_capacity[0] * sizeof( size_t ) );
    }
    if ( d_size[1] + N > d_capacity[1] ) {
        while ( d_size[1] + N > d_capacity[1] )
            d_capacity[1] = std::max<size_t>( 1024, 2 * d_capacity[1] );
        d_remote = (size_t *) realloc( d_remote, d_capacity[1] * sizeof( size_t ) );
    }
}


/****************************************************************
 * Get the number of non-zero columns in a row                   *
 ****************************************************************/
std::array<size_t, 2> GetRowHelper::NNZ( size_t row ) const
{
    const size_t begin = d_leftDOF->beginDOF();
    const size_t end   = d_leftDOF->endDOF();
    AMP_ASSERT( row >= begin && row < end );
    return d_NNZ[row - begin];
}


/****************************************************************
 * Get a row                                                     *
 ****************************************************************/
void GetRowHelper::getRow2( size_t row, size_t *local, size_t *remote ) const
{
    const size_t begin = d_leftDOF->beginDOF();
    const size_t end   = d_leftDOF->endDOF();
    AMP_ASSERT( row >= begin && row < end );
    auto N        = d_NNZ[row - begin];
    auto p_local  = &d_local[d_localOffset[row - begin]];
    auto p_remote = &d_remote[d_remoteOffset[row - begin]];
    if ( local ) {
        for ( size_t i = 0; i < N[0]; i++ )
            local[i] = p_local[i];
    }
    if ( remote ) {
        for ( size_t i = 0; i < N[1]; i++ )
            remote[i] = p_remote[i];
    }
}


} // namespace AMP::LinearAlgebra
