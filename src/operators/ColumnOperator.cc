#include "operators/ColumnOperator.h"
#include "utils/Utilities.h"
#include "vectors/MultiVariable.h"
#include "ProfilerApp.h"

namespace AMP {
  namespace Operator {

    void
      ColumnOperator :: apply(AMP::LinearAlgebra::Vector::const_shared_ptr f,
          AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r,
          const double a, const double b)
      {
        for(unsigned int i = 0; i < d_Operators.size(); i++)
        {
          d_Operators[i]->apply(f, u, r, a, b);
        }
      }

    AMP::shared_ptr<OperatorParameters>
      ColumnOperator :: getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u)
      {
        AMP::shared_ptr<AMP::Database> db;
        AMP::shared_ptr<ColumnOperatorParameters> opParameters(new ColumnOperatorParameters(db));

        (opParameters->d_OperatorParameters).resize(d_Operators.size());

        for(unsigned int i = 0; i < d_Operators.size(); i++)
        {
          (opParameters->d_OperatorParameters)[i] = (d_Operators[i]->getJacobianParameters(u));
        }
        return opParameters;
      }

    void
      ColumnOperator :: reset(const AMP::shared_ptr<OperatorParameters>& params)
      {
        AMP::shared_ptr<ColumnOperatorParameters> columnParameters =
          AMP::dynamic_pointer_cast<ColumnOperatorParameters>(params);

        AMP_INSIST( (columnParameters.get() != NULL), "ColumnOperator::reset parameter object is NULL" );

        AMP_INSIST( ( ((columnParameters->d_OperatorParameters).size()) == (d_Operators.size()) ), " std::vector sizes do not match! " );

        for(unsigned int i = 0; i < d_Operators.size(); i++)
        {
          d_Operators[i]->reset((columnParameters->d_OperatorParameters)[i]);
        }
      }

    void
      ColumnOperator :: append(AMP::shared_ptr< Operator > op)
      {
        AMP_INSIST( (op.get() != NULL), "AMP::ColumnOperator::appendRow input argument is a NULL operator");

        d_Operators.push_back(op);
      }

    AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getInputVariable()
    {
      AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable( new AMP::LinearAlgebra::MultiVariable("ColumnVariable"));

      for(unsigned int i = 0; i < d_Operators.size(); i++)
      {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = d_Operators[i]->getInputVariable();
        if(opVar.get()!=NULL)
        {
          retVariable->add(opVar);
        }
      }
      retVariable->removeDuplicateVariables();

      return retVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr ColumnOperator::getOutputVariable()
    {
      AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable( new AMP::LinearAlgebra::MultiVariable("ColumnVariable"));

      for(unsigned int i = 0; i < d_Operators.size(); i++)
      {
        AMP::LinearAlgebra::Variable::shared_ptr opVar = d_Operators[i]->getOutputVariable();
        if(opVar.get()!=NULL)
        {
          retVariable->add(opVar);
        }
      }
      retVariable->removeDuplicateVariables();

      return retVariable;
    }

    bool
      ColumnOperator::isValidInput(AMP::shared_ptr<AMP::LinearAlgebra::Vector> &u)
      {
        bool bRetVal=true;

        for(unsigned int i = 0; i < d_Operators.size(); i++)
        {
          bRetVal = bRetVal && d_Operators[i]->isValidInput(u);
        }

        return bRetVal;
      }

  }
}

