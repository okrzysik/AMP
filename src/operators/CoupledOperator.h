#ifndef included_AMP_CoupledOperator
#define included_AMP_CoupledOperator

#include "ColumnOperator.h"
#include "CoupledOperatorParameters.h"
#include "vectors/Vector.h"
#include "utils/Utilities.h"
#include <vector>

namespace AMP {
namespace Operator {


/**
   A class for representing a coupled operator combining a NodeToGaussPointOperator,
   a CopyOperator, a MapOperator, and a d_BVPOperator.
  */
class CoupledOperator : public ColumnOperator
{
public :
    CoupledOperator(const AMP::shared_ptr<OperatorParameters>& params);

    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
			AMP::LinearAlgebra::Vector::shared_ptr f ) override;


    /**
     * Default base class implementation of the residual: f-L(u)
     * \param f: shared pointer to const vector rhs
     * \param u: shared pointer to const vector u
     * \param r: shared pointer to vector residual
     */
    virtual void residual(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
			  AMP::LinearAlgebra::Vector::const_shared_ptr u, 
			  AMP::LinearAlgebra::Vector::shared_ptr r) override;


    AMP::shared_ptr<AMP::Operator::Operator> getMapOperator() {
        return d_Operators[2];
    }

    void setMapOperator( AMP::shared_ptr<AMP::Operator::Operator> op ) {
        d_Operators[2] = op;
    }

    AMP::shared_ptr<AMP::Operator::Operator> getBVPOperator() {
        return d_Operators[3];
    }

    void setBVPOperator( AMP::shared_ptr<AMP::Operator::Operator> op ) {
        d_Operators[3] = op;
    }

    virtual AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() {
        return d_Operators[3]->getOutputVariable();
    }

    virtual void append(AMP::shared_ptr< Operator > op) {
        AMP_ASSERT(d_Operators.size() < 4);
        AMP_ASSERT(op.get() != NULL);
        d_Operators.push_back(op);
    }

    bool isValidInput(AMP::shared_ptr<AMP::LinearAlgebra::Vector> &u)
    {
        return d_Operators[3]->isValidInput(u);
    }

    void setFrozenGaussPointVector(AMP::LinearAlgebra::Vector::shared_ptr u) {
        d_frozenGaussPointVector = u;
    }

    AMP::shared_ptr<OperatorParameters> getParameters(const std::string type,
                     AMP::LinearAlgebra::Vector::const_shared_ptr u,
                     AMP::shared_ptr<OperatorParameters> params = NULL ) override
    {
       return (d_Operators[3]->getParameters(type, u, params));
    }

    virtual ~CoupledOperator() { }

protected:

    AMP::LinearAlgebra::Vector::shared_ptr d_frozenGaussPointVector;

};


}
}

#endif

