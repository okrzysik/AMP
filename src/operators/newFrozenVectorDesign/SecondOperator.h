#ifndef included_AMP_SecondOperator
#define included_AMP_SecondOperator

#include "AMP/operators/newFrozenVectorDesign/OnePointOperator.h"
#include "AMP/vectors/MultiVariable.h"

namespace AMP {
namespace Operator {

class SecondOperator : public OnePointOperator
{
public:
    explicit SecondOperator( const std::shared_ptr<OperatorParameters> &params )
        : OnePointOperator( params )
    {
        d_constant = 3.0;
        d_primaryVar.reset(
            new AMP::LinearAlgebra::Variable( params->d_db->getString( "PrimaryVariable" ) ) );
        d_secondaryVar.reset(
            new AMP::LinearAlgebra::Variable( params->d_db->getString( "SecondaryVariable" ) ) );
    }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override
    {
        AMP::LinearAlgebra::Vector::const_shared_ptr inP =
            u->constSubsetVectorForVariable( d_primaryVar );
        AMP::LinearAlgebra::Vector::const_shared_ptr inS =
            u->constSubsetVectorForVariable( d_secondaryVar );
        AMP::LinearAlgebra::Vector::shared_ptr out = r->subsetVectorForVariable( d_primaryVar );
        out->linearSum( d_constant, *inP, 1.0, *inS );
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override
    {
        std::shared_ptr<AMP::LinearAlgebra::MultiVariable> retVariable(
            new AMP::LinearAlgebra::MultiVariable( "MultiVariable" ) );
        retVariable->add( d_primaryVar );
        retVariable->add( d_secondaryVar );
        return retVariable;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override { return d_primaryVar; }

protected:
    AMP::LinearAlgebra::Variable::shared_ptr d_primaryVar;
    AMP::LinearAlgebra::Variable::shared_ptr d_secondaryVar;

private:
};
} // namespace Operator
} // namespace AMP

#endif
