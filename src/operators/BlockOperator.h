
#ifndef included_AMP_BlockOperator
#define included_AMP_BlockOperator

#include "AMP/operators/Operator.h"

#include <vector>

namespace AMP::Operator {

class BlockOperator : public Operator
{
public:
    BlockOperator();

    explicit BlockOperator( std::shared_ptr<const OperatorParameters> params );

    virtual ~BlockOperator() {}

    std::string type() const override { return "BlockOperator"; }

    void setNumRowBlocks( int val );

    void setNumColumnBlocks( int val );

    void allocateBlocks();

    bool supportsMatrixFunctions();

    void setBlock( int row, int col, std::shared_ptr<Operator> op );

    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable() const override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getInputVariable() const override;

    void computeFirstIndices();

    int getNumRows();

    int getNumColumns();

    int getNumRowsForBlock( int id );

    int getNumColumnsForBlock( int id );

    static void
    getRow( void *object, int row, std::vector<size_t> &cols, std::vector<double> &values );

    void getRowForBlock( int locRow,
                         int blkRowId,
                         int blkColId,
                         std::vector<size_t> &locCols,
                         std::vector<double> &values );

protected:
    int d_iNumRowBlocks;
    int d_iNumColumnBlocks;

    std::vector<std::vector<std::shared_ptr<Operator>>> d_blocks;

    std::vector<int> d_firstRowId;
    std::vector<int> d_firstColumnId;

private:
};
} // namespace AMP::Operator

#endif
