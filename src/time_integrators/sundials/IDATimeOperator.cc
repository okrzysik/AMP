#include "time_integrators/sundials/IDATimeOperator.h"


namespace AMP{
namespace TimeIntegrator{
    
    IDATimeOperator::IDATimeOperator(boost::shared_ptr<AMP::Operator::OperatorParameters > in_params):TimeOperator(in_params)
    {
        
      d_cloningHappened = false;
      // BP, commenting out because this class does not define
      // a reset and the reset called is the reset for the base
      // class which already has been called once
      //        reset(in_params);
      
      //JL
      d_beta = 1.0;
      
    }
    
    IDATimeOperator::~IDATimeOperator()
    {
    }


    /*
    void
    IDATimeOperator::reset(const boost::shared_ptr<OperatorParameters>& in_params)
    {
        boost::shared_ptr<IDATimeOperatorParameters> params = boost::dynamic_pointer_cast<IDATimeOperatorParameters>(in_params);
        
        getFromInput(params->d_db);
        
    }
    */

    void
    IDATimeOperator::apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
                           AMP::LinearAlgebra::Vector::const_shared_ptr u,
                           AMP::LinearAlgebra::Vector::shared_ptr r,
                           const double a, const double b)
    {
        if (d_cloningHappened==0)
        {
          d_pScratchVector = r->cloneVector();
          d_pScratchVector->zero();
          d_cloningHappened = 1;
        }
    
        d_pMassOperator->apply(f, d_pIDATimeDerivative, d_pScratchVector, 1.0, 0.0);

        if(d_iDebugPrintInfoLevel>4)
          {
            AMP::pout << "Output of M * yp in IDATimeOperator" << std::endl;
            AMP::pout << d_pScratchVector << std::endl;
          }
        
        if(d_pAlgebraicVariable.get()!=NULL)          
          {
            boost::shared_ptr<AMP::LinearAlgebra::Vector> algebraicComponent = d_pScratchVector->subsetVectorForVariable(d_pAlgebraicVariable);
            algebraicComponent->zero();
          }
        
        d_pRhsOperator->apply(d_pScratchVector, u, r, 1.0, 1.0);

        bool dpSourceTermNull = (d_pSourceTerm.get()==NULL);
    
        if(d_iDebugPrintInfoLevel>5)
          {
            AMP::pout << "Output of M * yp-frhs(y,t) in IDATimeOperator" << std::endl;
            AMP::pout << r << std::endl;
          }
        
        if(!dpSourceTermNull)
          {
            r->axpby(-1.0, 1.0, d_pSourceTerm);
          }
        if(f.get()!=NULL)
          {
            r->axpby(b, a, f);
          }
        else
          {
            r->scale(a);
          }

        if(d_iDebugPrintInfoLevel>6)
          {
            AMP::pout << "Output of M * yp-frhs(y,t)-g in IDATimeOperator" << std::endl;
            AMP::pout << r << std::endl;
          }
        
    }

}
}

