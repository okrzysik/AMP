
#include "NonlinearBVPOperator.h"

#include "utils/Utilities.h"

#include <stdexcept>

namespace AMP {
  namespace Operator {

    NonlinearBVPOperator :: NonlinearBVPOperator(const boost::shared_ptr<BVPOperatorParameters>& params)
      : Operator (params) 
    {
      d_volumeOperator = params->d_volumeOperator;
      d_boundaryOperator = params->d_boundaryOperator;
      d_Mesh = d_volumeOperator->getMesh();
    }

    void
      NonlinearBVPOperator :: apply(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &f, 
          const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector>&r, double a, double b)
      {
        boost::shared_ptr<AMP::LinearAlgebra::Vector> nullRhs;

        AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );

        AMP::LinearAlgebra::Variable::shared_ptr outputVariable = this->getOutputVariable();
        AMP::LinearAlgebra::Vector::shared_ptr rInternal = r->subsetVectorForVariable( outputVariable );

        AMP_INSIST( ((rInternal.get()) != NULL), "NULL Internal Residual Vector" );

        d_volumeOperator->apply(nullRhs, u, r, 1.0, 0.0);

        if(d_iDebugPrintInfoLevel>3)
        {
          AMP::pout << "L2 Norm of r in NonlinearBVPOperator volumeOperator::apply is : "
            << r->L2Norm() << std::endl;
        }

        if(d_iDebugPrintInfoLevel>3)
        {
          AMP::pout << "L2 Norm of rInternal in NonlinearBVPOperator volumeOperator::apply is : "
            << rInternal->L2Norm() << std::endl;
        }

        d_boundaryOperator->apply(nullRhs, u, r, 1.0, 0.0);

        if(d_iDebugPrintInfoLevel>3)
        {
          AMP::pout << "L2 Norm of r in NonlinearBVPOperator boundaryOperator::apply is : "
            << r->L2Norm() << std::endl;
        }

        if(d_iDebugPrintInfoLevel>3)
        {
          AMP::pout << "L2 Norm of rInternal in NonlinearBVPOperator boundaryOperator::apply is : "
            << rInternal->L2Norm() << std::endl;
        }

        if(f.get() == NULL) {
          rInternal->scale(a);
        } else {
          AMP::LinearAlgebra::Vector::shared_ptr fInternal = f->subsetVectorForVariable( this->getOutputVariable() );
          if(fInternal.get() == NULL) {
            rInternal->scale(a);
          } else {
            rInternal->axpby(b, a, fInternal);
          }
        }

        if(d_iDebugPrintInfoLevel>2)
        {
          AMP::pout << "L2 Norm of output of NonlinearBVPOperator :: apply is : " 
            << rInternal->L2Norm() << std::endl;
        }
      }

    void NonlinearBVPOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      boost::shared_ptr<BVPOperatorParameters> inParams = 
        boost::dynamic_pointer_cast<BVPOperatorParameters>(params);

      AMP_INSIST( (inParams.get() != NULL), "NonlinearBVPOperator :: reset Null parameter" );

      d_volumeOperator->reset(inParams->d_volumeOperatorParams);
      d_boundaryOperator->reset(inParams->d_boundaryOperatorParams);
    }

    boost::shared_ptr<OperatorParameters> NonlinearBVPOperator :: getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u)
    {
      boost::shared_ptr<AMP::Database> db;
      boost::shared_ptr<BVPOperatorParameters> outParams(new BVPOperatorParameters(db));

      outParams->d_volumeOperatorParams = d_volumeOperator->getJacobianParameters(u);
      outParams->d_boundaryOperatorParams = d_boundaryOperator->getJacobianParameters(u);

      return outParams;
    }

    void NonlinearBVPOperator :: modifyRHSvector(AMP::LinearAlgebra::Vector::shared_ptr rhs) {
      (this->getBoundaryOperator())->addRHScorrection(rhs);
      (this->getBoundaryOperator())->setRHScorrection(rhs);
    }

    void NonlinearBVPOperator :: modifyInitialSolutionVector(AMP::LinearAlgebra::Vector::shared_ptr sol) {
      (this->getBoundaryOperator())->modifyInitialSolutionVector(sol);
    }

  }
}


