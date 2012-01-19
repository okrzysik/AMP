
#ifndef included_AMP_BlockOperator
#define included_AMP_BlockOperator

#include "operators/Operator.h"

#include <vector>

namespace AMP {
  namespace Operator {

    class BlockOperator : public Operator 
    {
      public:

        BlockOperator();

        BlockOperator(const boost::shared_ptr<OperatorParameters>& params);

        ~BlockOperator() { }

        void setNumRowBlocks(int val);

        void setNumColumnBlocks(int val);

        void allocateBlocks();

        bool supportsMatrixFunctions();

        void setBlock(int row, int col, boost::shared_ptr<Operator> op);

        void reset(const boost::shared_ptr<OperatorParameters>& params);

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

        void computeFirstIndices();

        int getNumRows();

        int getNumColumns();

        int getNumRowsForBlock(int id);

        int getNumColumnsForBlock(int id);

        static void getRow(void * object, int row, std::vector<unsigned int> &cols, std::vector<double> &values);

        void getRowForBlock(int locRow, int blkRowId, int blkColId,
            std::vector<unsigned int> &locCols, std::vector<double> &values);

      protected :

        int d_iNumRowBlocks;
        int d_iNumColumnBlocks;

        std::vector<std::vector<boost::shared_ptr<Operator> > > d_blocks;

        std::vector<int> d_firstRowId;
        std::vector<int> d_firstColumnId;

      private :

    };

  }
}

#endif


