#ifndef included_AMP_RowOperator
#define included_AMP_RowOperator


#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/vectors/Vector.h"

#include <memory>


namespace AMP::Operator {

class RowOperator : public Operator
{
public:
    explicit RowOperator( std::shared_ptr<const OperatorParameters> params );

    virtual ~RowOperator() = default;

    std::string type() const override { return "RowOperator"; }

    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    void resetScaling( int idx, double a ) { scalea[idx] = a; }

    void append( std::shared_ptr<Operator> op, double a );

    void setJacobianFlag() { getAllJacobian = true; }

    void setJacobianParametersSize( const int paramSz ) { d_paramsize = paramSz; }
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    std::shared_ptr<AMP::LinearAlgebra::Variable> getOutputVariable() const override
    {
        return d_Operators[0]->getOutputVariable();
    }

    std::shared_ptr<Operator> getOperator( const int i ) { return d_Operators[i]; }

    int getNumberOfOperators( void ) { return d_Operators.size(); }

protected:
    std::shared_ptr<OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   std::shared_ptr<OperatorParameters> params = nullptr ) override;


    std::vector<std::shared_ptr<Operator>> d_Operators;

    std::vector<double> scalea;
    int d_paramsize;

    bool getAllJacobian;

private:
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_OutputVariable;
};


} // namespace AMP::Operator

#endif
