
#ifndef included_PetscMatrixShellOperator
#define included_PetscMatrixShellOperator

#include "LinearOperator.h"
#include "utils/AMP_MPI.h"

extern "C" {

#ifdef MPICH_SKIP_MPICXX
#define _FIX_FOR_PETSC_MPI_CXX
#undef MPICH_SKIP_MPICXX
#endif

#ifdef OMPI_SKIP_MPICXX
#define _FIX_FOR_PETSC_OMPI_CXX
#undef OMPI_SKIP_MPICXX
#endif

#include "petscmat.h"

#ifdef _FIX_FOR_PETSC_MPI_CXX
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif
#endif

#ifdef _FIX_FOR_PETSC_OMPI_CXX
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#endif

}


namespace AMP {
  namespace Operator {

    class PetscMatrixShellOperator : public LinearOperator {
      public:

        PetscMatrixShellOperator(const boost::shared_ptr<OperatorParameters>& params);

        ~PetscMatrixShellOperator() { }

        static PetscErrorCode mult(Mat, Vec, Vec);

        void setMatLocalRowSize(int val); 

        void setMatLocalColumnSize(int val);  

        void setComm(AMP_MPI comm);

        void setOperator(boost::shared_ptr<Operator> op);       

        void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, AMP::LinearAlgebra::Vector::const_shared_ptr u,
            AMP::LinearAlgebra::Vector::shared_ptr r, const double a = -1.0, const double b = 1.0);

        void reset(const boost::shared_ptr<OperatorParameters>& params);

        AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable();

        AMP::LinearAlgebra::Variable::shared_ptr getInputVariable();

      private:

        boost::shared_ptr<Operator> d_operator;
        int d_iMatLocalRowSize;
        int d_iMatLocalColumnSize;
        Mat d_mat;
        AMP_MPI d_comm;
    };

  }
}

#endif


