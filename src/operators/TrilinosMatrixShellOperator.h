
#ifndef included_TrilinosMatrixShellOperator
#define included_TrilinosMatrixShellOperator

#include "LinearOperator.h"

#include "ampmesh/MeshManager.h"

#include "ml_include.h"

namespace AMP {
  namespace Operator {

    class TrilinosMatrixShellOperator : public LinearOperator {
      public:

        TrilinosMatrixShellOperator(const boost::shared_ptr<OperatorParameters>& params);

        ~TrilinosMatrixShellOperator() { }

        void setMeshManager(AMP::Mesh::MeshManager::shared_ptr manager);

        void setOperator(boost::shared_ptr<Operator> op);     

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
            AMP::LinearAlgebra::Vector::shared_ptr &r, const double a = -1.0, const double b = 1.0);

        void reset(const boost::shared_ptr<OperatorParameters>& params);

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable(int varId = -1);

        static int matVec(ML_Operator *data, int in_length, double in[], int out_length, double out[]);

        static int getRow(ML_Operator *data, int N_requested_rows, int requested_rows[],
            int allocated_space, int columns[], double values[], int row_lengths[] );

        void setGetRow(void (*func)(void* object, int row, std::vector<unsigned int> &cols, std::vector<double> &values));

        void getColumn(int column, std::vector<unsigned int> &rows, std::vector<double> &values);

        size_t getMatrixSize();

      private:

        AMP::Mesh::MeshManager::shared_ptr d_meshManager;

        boost::shared_ptr<Operator> d_operator;

        void (*d_getRow)(void* object, int row, std::vector<unsigned int> &cols, std::vector<double> &values);

    };

  }
}

#endif


