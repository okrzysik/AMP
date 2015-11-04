#ifndef included_AMP_CoupledFlowFrapconOperator
#define included_AMP_CoupledFlowFrapconOperator

#include "operators/ColumnOperator.h"
#include "operators/subchannel/CoupledFlowFrapconOperatorParameters.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include <vector>

namespace AMP {
namespace Operator {

  class CoupledFlowFrapconOperator : public ColumnOperator
  {
    public :
      CoupledFlowFrapconOperator(const AMP::shared_ptr<OperatorParameters>& params);

      AMP::shared_ptr<OperatorParameters>
        getJacobianParameters(const AMP::LinearAlgebra::Vector::shared_ptr & u)
        {
          return (d_Operators[2]->getJacobianParameters(u));
        }

      void reset(const AMP::shared_ptr<OperatorParameters>& params){
        d_Operators[2]->reset(params);
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() {
        return d_Operators[4]->getOutputVariable();
      }

      virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_Operators[4]->getOutputVariable();
      }

      void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
		  AMP::LinearAlgebra::Vector::shared_ptr f) override;

      virtual void append(AMP::shared_ptr< Operator > op) {
        AMP_ASSERT(d_Operators.size() < 3);
        AMP_ASSERT(op.get() != NULL);
        d_Operators.push_back(op);
      }

      virtual ~CoupledFlowFrapconOperator() { }

    protected :

    private :

      AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

      int d_numpoints; /**< Number of points in z direction */

      std::vector<double> d_zPoints; /**< vector to hold z locations */
      AMP::LinearAlgebra::Vector::shared_ptr d_flowInput; 
      AMP::LinearAlgebra::Vector::shared_ptr d_flowOutput; 
      AMP::LinearAlgebra::Vector::shared_ptr d_localCladVec;

      AMP::Operator::MapOperator::shared_ptr d_flowInternal3to1;
      AMP::Operator::MapOperator::shared_ptr d_flowInternal1to3;
  };

}
}

#endif

