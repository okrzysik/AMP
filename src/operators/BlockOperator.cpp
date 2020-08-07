#include "AMP/operators/BlockOperator.h"
#include "AMP/operators/BlockOperatorParameters.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/vectors/MultiVariable.h"

#include <algorithm>

namespace AMP {
namespace Operator {

BlockOperator::BlockOperator() : Operator()
{
    d_iNumRowBlocks    = -1234;
    d_iNumColumnBlocks = -5678;
}

BlockOperator::BlockOperator( const std::shared_ptr<OperatorParameters> & ) : Operator()
{
    d_iNumRowBlocks    = -1234;
    d_iNumColumnBlocks = -5678;
}

void BlockOperator::reset( const std::shared_ptr<OperatorParameters> &params )
{
    std::shared_ptr<BlockOperatorParameters> myParams =
        std::dynamic_pointer_cast<BlockOperatorParameters>( params );
    for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
        for ( int j = 0; j < d_iNumColumnBlocks; j++ ) {
            d_blocks[i][j]->reset( ( myParams->d_blockParams )[i][j] );
        }
    }
}

AMP::LinearAlgebra::Variable::shared_ptr BlockOperator::getOutputVariable()
{
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> var(
        new AMP::LinearAlgebra::MultiVariable( "BlockVariable" ) );
    for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
        var->add( d_blocks[i][0]->getOutputVariable() );
    }
    var->removeDuplicateVariables();

    return var;
}

AMP::LinearAlgebra::Variable::shared_ptr BlockOperator::getInputVariable()
{
    std::shared_ptr<AMP::LinearAlgebra::MultiVariable> var(
        new AMP::LinearAlgebra::MultiVariable( "BlockVariable" ) );
    for ( int i = 0; i < d_iNumColumnBlocks; i++ ) {
        var->add( d_blocks[0][i]->getInputVariable() );
    }
    var->removeDuplicateVariables();

    return var;
}

void BlockOperator::allocateBlocks()
{
    d_blocks.resize( d_iNumRowBlocks );

    for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
        d_blocks[i].resize( d_iNumColumnBlocks );
    }
}

void BlockOperator::setBlock( int row, int col, std::shared_ptr<Operator> op )
{
    d_blocks[row][col] = op;
}

void BlockOperator::setNumRowBlocks( int val ) { d_iNumRowBlocks = val; }

void BlockOperator::setNumColumnBlocks( int val ) { d_iNumColumnBlocks = val; }

bool BlockOperator::supportsMatrixFunctions()
{
    for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
        for ( int j = 0; j < d_iNumColumnBlocks; j++ ) {
            std::shared_ptr<BlockOperator> blockOp =
                std::dynamic_pointer_cast<BlockOperator>( d_blocks[i][j] );
            if ( blockOp == nullptr ) {
                std::shared_ptr<LinearOperator> matOp =
                    std::dynamic_pointer_cast<LinearOperator>( d_blocks[i][j] );
                if ( matOp == nullptr ) {
                    return false;
                }
            } else {
                if ( !( blockOp->supportsMatrixFunctions() ) ) {
                    return false;
                }
            }
        }
    }
    return true;
}

void BlockOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                           AMP::LinearAlgebra::Vector::shared_ptr r )
{

    AMP::LinearAlgebra::Variable::shared_ptr tmpOutVar = getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable( tmpOutVar );

    rInternal->zero(rInternal);
    AMP::LinearAlgebra::Vector::shared_ptr rCopy = rInternal->cloneVector();
    for ( int j = 0; j < d_iNumColumnBlocks; j++ ) {
        for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
            AMP::LinearAlgebra::Vector::shared_ptr nullVec;
            d_blocks[i][j]->apply( u, rCopy );
        }
        rInternal->add( rInternal, rCopy, rInternal );
    }
}

void BlockOperator::computeFirstIndices()
{
    d_firstRowId.resize( d_iNumRowBlocks, 0 );
    for ( int i = 1; i < d_iNumRowBlocks; i++ ) {
        d_firstRowId[i] = d_firstRowId[i - 1] + getNumRowsForBlock( i - 1 );
    }

    d_firstColumnId.resize( d_iNumColumnBlocks, 0 );
    for ( int i = 1; i < d_iNumColumnBlocks; i++ ) {
        d_firstColumnId[i] = d_firstColumnId[i - 1] + getNumColumnsForBlock( i - 1 );
    }
}

int BlockOperator::getNumRows()
{
    int result = 0;
    for ( int i = 0; i < d_iNumRowBlocks; i++ ) {
        result += getNumRowsForBlock( i );
    }
    return result;
}

int BlockOperator::getNumColumns()
{
    int result = 0;
    for ( int i = 0; i < d_iNumColumnBlocks; i++ ) {
        result += getNumColumnsForBlock( i );
    }
    return result;
}

int BlockOperator::getNumRowsForBlock( int id )
{
    int result = 0;
    std::shared_ptr<BlockOperator> blockOp =
        std::dynamic_pointer_cast<BlockOperator>( d_blocks[id][0] );
    if ( blockOp == nullptr ) {
        std::shared_ptr<LinearOperator> matOp =
            std::dynamic_pointer_cast<LinearOperator>( d_blocks[id][0] );
        if ( matOp == nullptr ) {
            AMP_ERROR( "This is not supported." );
        } else {
            result = ( ( matOp->getMatrix() )->numGlobalRows() );
        }
    } else {
        result = ( blockOp->getNumRows() );
    }
    return result;
}

int BlockOperator::getNumColumnsForBlock( int id )
{
    int result = 0;
    std::shared_ptr<BlockOperator> blockOp =
        std::dynamic_pointer_cast<BlockOperator>( d_blocks[0][id] );
    if ( blockOp == nullptr ) {
        std::shared_ptr<LinearOperator> matOp =
            std::dynamic_pointer_cast<LinearOperator>( d_blocks[0][id] );
        if ( matOp == nullptr ) {
            AMP_ERROR( "This is not supported." );
        } else {
            result = ( ( matOp->getMatrix() )->numGlobalColumns() );
        }
    } else {
        result = ( blockOp->getNumColumns() );
    }
    return result;
}

void BlockOperator::getRow( void *object,
                            int row,
                            std::vector<size_t> &cols,
                            std::vector<double> &values )
{
    auto *myObject = static_cast<BlockOperator *>( object );

    AMP_ASSERT( row >= 0 );
    AMP_ASSERT( row < ( myObject->getNumRows() ) );

    auto blkPtr =
        std::lower_bound( myObject->d_firstRowId.begin(), myObject->d_firstRowId.end(), row );

    int blkRowId = -1;
    if ( blkPtr == myObject->d_firstRowId.end() ) {
        blkRowId = myObject->d_iNumRowBlocks - 1;
    } else if ( ( *blkPtr ) > row ) {
        blkRowId = blkPtr - myObject->d_firstRowId.begin() - 1;
    } else {
        AMP_ASSERT( ( *blkPtr ) == row );
        blkRowId = blkPtr - myObject->d_firstRowId.begin();
    }
    AMP_ASSERT( blkRowId >= 0 );
    AMP_ASSERT( blkRowId < myObject->d_iNumRowBlocks );
    AMP_ASSERT( myObject->d_firstRowId[blkRowId] <= row );
    if ( blkRowId < ( myObject->d_iNumRowBlocks - 1 ) ) {
        AMP_ASSERT( myObject->d_firstRowId[blkRowId + 1] > row );
    }

    int locRow = row - myObject->d_firstRowId[blkRowId];

    cols.clear();
    values.clear();
    for ( int blkColId = 0; blkColId < myObject->d_iNumColumnBlocks; blkColId++ ) {
        std::vector<size_t> tmpCols;
        std::vector<double> tmpVals;
        myObject->getRowForBlock( locRow, blkRowId, blkColId, tmpCols, tmpVals );
        for ( auto &tmpCol : tmpCols ) {
            tmpCol += myObject->d_firstColumnId[blkColId];
        }
        cols.insert( cols.end(), tmpCols.begin(), tmpCols.end() );
        values.insert( values.end(), tmpVals.begin(), tmpVals.end() );
    }
}

void BlockOperator::getRowForBlock( int locRow,
                                    int blkRowId,
                                    int blkColId,
                                    std::vector<size_t> &locCols,
                                    std::vector<double> &values )
{
    std::shared_ptr<BlockOperator> blockOp =
        std::dynamic_pointer_cast<BlockOperator>( d_blocks[blkRowId][blkColId] );
    if ( blockOp == nullptr ) {
        std::shared_ptr<LinearOperator> matOp =
            std::dynamic_pointer_cast<LinearOperator>( d_blocks[blkRowId][blkColId] );
        if ( matOp == nullptr ) {
            AMP_ERROR( "This is not supported." );
        } else {
            AMP::LinearAlgebra::Matrix::shared_ptr mat = matOp->getMatrix();
            if ( mat == nullptr ) {
                AMP_ERROR( "Matrix is NULL" );
            } else {
                mat->getRowByGlobalID( locRow, locCols, values );
            }
        }
    } else {
        getRow( blockOp.get(), locRow, locCols, values );
    }
}
} // namespace Operator
} // namespace AMP
