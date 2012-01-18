
#ifndef included_AMP_DirichletMatrixCorrection
#define included_AMP_DirichletMatrixCorrection

#include "BoundaryOperator.h"
#include "DirichletMatrixCorrectionParameters.h"
#include "DirichletVectorCorrection.h"

#include "vectors/Variable.h"

namespace AMP {
  namespace Operator {

    /**
      A class used to impose Dirichlet boundary conditions for a linear operator. For a linear operator, imposing
      Dirichlet boundary conditions involves the following steps:
      1) Modify the entries of the matrix appropriately.
      2) Add a vector of corrections to the RHS vector
      3) Set the Dirichlet values at the appropriate locations in the RHS vector.
      */
    class DirichletMatrixCorrection : public BoundaryOperator
    {
      public :

        /**
          Constructor
          */
        DirichletMatrixCorrection(const boost::shared_ptr<DirichletMatrixCorrectionParameters> & params)
          : BoundaryOperator (params)
        {
          d_variable = params->d_variable;
          reset(params);
        }

        /**
          Destructor
          */
        ~DirichletMatrixCorrection() { }

        /**
          Set the variable for the vector that will used with this operator.
          */
        void setVariable(const AMP::LinearAlgebra::Variable::shared_ptr & var) {
          d_variable = var;
        }

        void apply(const AMP::LinearAlgebra::Vector::shared_ptr &,
            const AMP::LinearAlgebra::Vector::shared_ptr &, AMP::LinearAlgebra::Vector::shared_ptr &,
            const double, const double)
        {
          //Do Nothing
        }

        void parseParams(const boost::shared_ptr<DirichletMatrixCorrectionParameters> & );

        void computeRHScorrection(const boost::shared_ptr<DirichletMatrixCorrectionParameters> & );

        /**
          This function modifies the entries of the matrix formed by the volume operator
          in order to impose Dirichlet boundary conditions. This function can also be used
          to change the Dirichlet boundary conditions, if required.
          */
        void reset(const boost::shared_ptr<OperatorParameters>& );

        /**
          Adds a vector to the RHS vector. This is one of the steps for imposing Dirichlet boundary conditions.  
          This step can be skipped if the Dirichlet boundary conditons are homogeneous.
          */
        void addRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
          if(!d_skipRHSaddCorrection) {
            AMP::LinearAlgebra::Vector::shared_ptr myRhs = rhs->subsetVectorForVariable(d_variable);
            myRhs->add(myRhs, d_rhsCorrectionAdd);
          }
        }

        /**
          Sets the Dirichlet values at the appropriate locations in the RHS vector. This is one
          of the steps for imposing Dirichlet boundary conditions.  
          */
        void setRHScorrection(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
          if(!d_skipRHSsetCorrection) {
            AMP::LinearAlgebra::Vector::shared_ptr emptyVec;
            d_rhsCorrectionSet->apply(emptyVec, emptyVec, rhs, 1.0, 0.0);
          }
        }

        std::vector<short int> getBoundaryIds() {
          return d_boundaryIds;
        }

        std::vector<std::vector<unsigned int> > getDofIds() {
          return d_dofIds;
        }

      protected :

        //This must be a simple variable not a dual or multivariable
        AMP::LinearAlgebra::Variable::shared_ptr d_variable;

        std::vector<short int> d_boundaryIds;

        std::vector<std::vector<double> > d_dirichletValues;

        std::vector<std::vector<unsigned int> > d_dofIds;

        AMP::LinearAlgebra::Vector::shared_ptr d_rhsCorrectionAdd;

        AMP::LinearAlgebra::Vector::shared_ptr d_dispVals;

        boost::shared_ptr<DirichletVectorCorrection> d_rhsCorrectionSet;

        bool d_symmetricCorrection;

        bool d_zeroDirichletBlock;

        bool d_skipRHSaddCorrection;

        bool d_skipRHSsetCorrection;

      private :

    };

  }
}

#endif

