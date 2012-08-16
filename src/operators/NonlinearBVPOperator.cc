#include "NonlinearBVPOperator.h"

#include "utils/Utilities.h"
#include "utils/ProfilerApp.h"
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
      NonlinearBVPOperator :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr f, 
          AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r, double a, double b)
      {
        PROFILE_START("apply");
        boost::shared_ptr<AMP::LinearAlgebra::Vector> nullRhs;

        AMP_INSIST( ((r.get()) != NULL), "NULL Residual Vector" );

        AMP::LinearAlgebra::Vector::shared_ptr rInternal = this->subsetOutputVector(r);

        AMP_INSIST( ((rInternal.get()) != NULL), "NULL Internal Residual Vector" );

        if(d_iDebugPrintInfoLevel>3)
        {
          AMP::pout << "L2 Norm of u in NonlinearBVPOperator volumeOperator::apply is : "
            << u->L2Norm() << std::endl;
        }

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
          AMP::LinearAlgebra::Vector::const_shared_ptr fInternal = this->subsetOutputVector(f);
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
        PROFILE_STOP("apply");
      }

    void NonlinearBVPOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      PROFILE_START("reset");
      boost::shared_ptr<BVPOperatorParameters> inParams = 
        boost::dynamic_pointer_cast<BVPOperatorParameters>(params);

      AMP_INSIST( (inParams.get() != NULL), "NonlinearBVPOperator :: reset Null parameter" );

      d_volumeOperator->reset(inParams->d_volumeOperatorParams);
      d_boundaryOperator->reset(inParams->d_boundaryOperatorParams);
      PROFILE_STOP("reset");
    }

    boost::shared_ptr<OperatorParameters> NonlinearBVPOperator :: getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u)
    {
      PROFILE_START("getJacobianParameters");
      boost::shared_ptr<AMP::Database> db;
      boost::shared_ptr<BVPOperatorParameters> outParams(new BVPOperatorParameters(db));

      outParams->d_volumeOperatorParams = d_volumeOperator->getJacobianParameters(u);
      outParams->d_boundaryOperatorParams = d_boundaryOperator->getJacobianParameters(u);

      PROFILE_STOP("getJacobianParameters");
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


