
#ifndef included_AMP_ColumnBoundaryOperator
#define included_AMP_ColumnBoundaryOperator

#include "AMP/operators/ColumnOperatorParameters.h"
#include "AMP/vectors/Vector.h"
#include "BoundaryOperator.h"

#include <vector>

namespace AMP {
namespace Operator {

/**
  A class for representing a composite boundary operator, F=(F1, F2, F3, .., Fk),
  where each of F1,.., Fk are boundary operators themselves. The user is expected to have
  created and initialized the operators F1,.., Fk. This class can be used to impose a mix of
  boundary
  conditions for a volume operator over the boundary
  */

typedef ColumnOperatorParameters ColumnBoundaryOperatorParameters;

class ColumnBoundaryOperator : public BoundaryOperator
{

public:
    explicit ColumnBoundaryOperator( const std::shared_ptr<OperatorParameters> &params )
        : BoundaryOperator( params )
    {
    }

    virtual ~ColumnBoundaryOperator() {}

    std::string type() const override { return "ColumnBoundaryOperator"; }

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr r ) override;

    std::shared_ptr<OperatorParameters>
    getParameters( const std::string &type,
                   AMP::LinearAlgebra::Vector::const_shared_ptr u,
                   std::shared_ptr<OperatorParameters> params = nullptr ) override;

    void reset( std::shared_ptr<const OperatorParameters> params ) override;

    /**
     * \param op
     *            shared pointer to an operator to append to the existing column of operators
     */
    virtual void append( std::shared_ptr<BoundaryOperator> op );

    std::shared_ptr<BoundaryOperator> getBoundaryOperator( int i ) { return d_Operators[i]; }

    size_t numberOfBoundaryOperators() { return d_Operators.size(); }

    void addRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr ) override;

    void setRHScorrection( AMP::LinearAlgebra::Vector::shared_ptr ) override;

    void modifyInitialSolutionVector( AMP::LinearAlgebra::Vector::shared_ptr ) override;

protected:
    std::vector<std::shared_ptr<BoundaryOperator>> d_Operators;

private:
};
} // namespace Operator
} // namespace AMP

#endif
