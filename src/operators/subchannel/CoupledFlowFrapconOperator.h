#ifndef included_AMP_CoupledFlowFrapconOperator
#define included_AMP_CoupledFlowFrapconOperator

#include "operators/ColumnOperator.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/subchannel/CoupledFlowFrapconOperatorParameters.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "utils/Utilities.h"
#include "vectors/Vector.h"
#include <vector>

namespace AMP {
namespace Operator {

class CoupledFlowFrapconOperator : public ColumnOperator
{
public:
    explicit CoupledFlowFrapconOperator( const AMP::shared_ptr<OperatorParameters> &params );

    virtual void reset( const AMP::shared_ptr<OperatorParameters> &params ) override
    {
        d_Operators[2]->reset( params );
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        return d_Operators[4]->getOutputVariable();
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_Operators[4]->getOutputVariable();
    }

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual void append( AMP::shared_ptr<Operator> op ) override
    {
        AMP_ASSERT( d_Operators.size() < 3 );
        AMP_ASSERT( op.get() != NULL );
        d_Operators.push_back( op );
    }

    AMP::shared_ptr<OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   AMP::shared_ptr<OperatorParameters> params = NULL ) override
    {
        return ( d_Operators[2]->getParameters( type, u, params ) );
    }

    virtual ~CoupledFlowFrapconOperator() {}

protected:
private:
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> d_SimpleVariable; /**< Simple Input Variable */

    int d_numpoints; /**< Number of points in z direction */

    std::vector<double> d_zPoints; /**< vector to hold z locations */
    AMP::LinearAlgebra::Vector::shared_ptr d_flowInput;
    AMP::LinearAlgebra::Vector::shared_ptr d_flowOutput;
    AMP::LinearAlgebra::Vector::shared_ptr d_localCladVec;

    AMP::Operator::MapOperator::shared_ptr d_flowInternal3to1;
    AMP::Operator::MapOperator::shared_ptr d_flowInternal1to3;
};
} // namespace Operator
} // namespace AMP

#endif
