#ifndef included_AMP_RowOperator
#define included_AMP_RowOperator

/* AMP files */
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "vectors/Vector.h"

#include "ColumnOperatorParameters.h"
#include "ampmesh/MeshManager.h"
#include "ampmesh/DOFMap.h"
#include "utils/Utilities.h"
#include "boost/shared_ptr.hpp"

namespace AMP{
namespace Operator {

  class RowOperator : public Operator
  {
    public :

      typedef boost::shared_ptr <RowOperator>  shared_ptr;

      RowOperator(const boost::shared_ptr<OperatorParameters> & params)
        : Operator () { 
          (void) params; 
          getAllJacobian = false;
          d_paramsize  = 1;
        }

      virtual ~RowOperator(){ }

      virtual void reset(const boost::shared_ptr<OperatorParameters>& params)
      {
        boost::shared_ptr<ColumnOperatorParameters> fParams =
          boost::dynamic_pointer_cast<ColumnOperatorParameters>(params);

        AMP_INSIST( (fParams.get() != NULL), "RowOperator::reset parameter object is NULL" );

        AMP_INSIST( ( ((fParams->d_OperatorParameters).size()) == (d_Operators.size()) ), " std::vector sizes do not match! " );

        for(unsigned int i = 0; i < d_Operators.size(); i++)
        {
          d_Operators[i]->reset((fParams->d_OperatorParameters)[i]);
        }
      }

      void resetScaling(int idx, double a)
      {
        scalea[idx]=a;
      }

      void  append(boost::shared_ptr< Operator > op, double a)
      {
        AMP_INSIST( (op.get() != NULL), "AMP::RowOperator::appendRow input argument is a NULL operator");

        d_Operators.push_back(op);
        scalea.push_back(a);
      }

      void setJacobianFlag()
      {
        getAllJacobian = true;
      }

      void setJacobianParametersSize(const int paramSz)
      {
        d_paramsize = paramSz;
      }
      virtual void apply(const AMP::LinearAlgebra::Vector::shared_ptr &f,
          const  AMP::LinearAlgebra::Vector::shared_ptr &u, AMP::LinearAlgebra::Vector::shared_ptr  &r,
          const double a = -1.0, const double b = 1.0);

      virtual boost::shared_ptr<OperatorParameters>
        getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u)
        {
          boost::shared_ptr<AMP::Database> db;

          boost::shared_ptr<ColumnOperatorParameters> opParameters(new ColumnOperatorParameters(db));
          
          boost::shared_ptr<OperatorParameters> rtParameters(new OperatorParameters(db));

          if(getAllJacobian)
          {
            (opParameters->d_OperatorParameters).resize(d_Operators.size());

            for(unsigned int i = 0; i < d_Operators.size(); i++)
            {
              (opParameters->d_OperatorParameters)[i] = (d_Operators[i]->getJacobianParameters(u));
            }

            rtParameters = boost::dynamic_pointer_cast<OperatorParameters>(opParameters);

          }else{
            (opParameters->d_OperatorParameters).resize(d_paramsize);

            for(int i = 0; i < d_paramsize; i++)
            {
              (opParameters->d_OperatorParameters)[i] = (d_Operators[i]->getJacobianParameters(u));
            }

            rtParameters = boost::dynamic_pointer_cast<OperatorParameters>(opParameters);
            //rtParameters = (d_Operators[0]->getJacobianParameters(u));

          }
          return rtParameters;
        }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        d_OutputVariable = d_Operators[0]->getOutputVariable();
        return d_OutputVariable;
      }

      boost::shared_ptr< Operator > getOperator(const int i){return d_Operators[i]; }

      int getNumberOfOperators(void){return d_Operators.size(); }

    protected :

      std::vector< boost::shared_ptr< Operator > > d_Operators;

      std::vector< double > scalea;
      int d_paramsize ;

      bool getAllJacobian;

    private:

      AMP::LinearAlgebra::Variable::shared_ptr d_OutputVariable;

  };

}
}

#endif

