
#ifndef included_AMP_BlockOperator
#define included_AMP_BlockOperator

#include "operators/Operator.h"

#include <vector>

namespace AMP {
namespace Operator {

class BlockOperator : public Operator {
public:
    BlockOperator();

    explicit BlockOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual ~BlockOperator() {}

    void setNumRowBlocks( int val );

    void setNumColumnBlocks( int val );

    void allocateBlocks();

    bool supportsMatrixFunctions();

    void setBlock( int row, int col, AMP::shared_ptr<Operator> op );

    void reset( const AMP::shared_ptr<OperatorParameters> &params );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

    void computeFirstIndices();

    int getNumRows();

    int getNumColumns();

    int getNumRowsForBlock( int id );

    int getNumColumnsForBlock( int id );

    static void
    getRow( void *object, int row, std::vector<unsigned int> &cols, std::vector<double> &values );

    void getRowForBlock( int locRow,
                         int blkRowId,
                         int blkColId,
                         std::vector<unsigned int> &locCols,
                         std::vector<double> &values );

protected:
    int d_iNumRowBlocks;
    int d_iNumColumnBlocks;

    std::vector<std::vector<AMP::shared_ptr<Operator>>> d_blocks;

    std::vector<int> d_firstRowId;
    std::vector<int> d_firstColumnId;

private:
};
}
}

#endif
