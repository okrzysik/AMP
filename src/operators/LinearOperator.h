
#ifndef included_AMP_LinearOperator
#define included_AMP_LinearOperator

#include "AMP/matrices/Matrix.h"
#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/vectors/Vector.h"
#include <memory>

namespace AMP {
namespace Operator {


/**
 * An abstract base class for representing a linear operator. This class
 * stores the matrix representation of the linear operator. It provides
 * an implementation of the apply() function.
 * @see Operator
 */
class LinearOperator : public Operator
{
public:
    /**
     * Constructor. This resets the matrix shared pointer.
     * @param [in] params
     */
    explicit LinearOperator( const std::shared_ptr<OperatorParameters> &params );

    //! Destructor
    virtual ~LinearOperator() {}

    /**
     * The apply function for this operator, A, performs the following operation:
     * f = A(u)
     * Here, A(u) is simply a Matrix-Vector multiplication.
     * @param [in] u input vector.
     * @param [out] f residual/output vector.
     */
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                        AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
     * @return The matrix representation of this linear operator.
     */
    std::shared_ptr<AMP::LinearAlgebra::Matrix> getMatrix();

    /**
     * Copies the shared pointer for the matrix representation of this linear operator.
     *  @param [in] in_mat The matrix representation of this linear operator.
     */
    virtual void setMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> in_mat );

    virtual void setVariables( AMP::LinearAlgebra::Variable::shared_ptr in,
                               AMP::LinearAlgebra::Variable::shared_ptr out )
    {
        d_inputVariable  = in;
        d_outputVariable = out;
    }

    AMP::LinearAlgebra::Variable::shared_ptr getInputVariable() override { return d_inputVariable; }
    AMP::LinearAlgebra::Variable::shared_ptr getOutputVariable() override
    {
        return d_outputVariable;
    }

protected:
    //! Empty constructor
    LinearOperator();

    // input and output variables
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_inputVariable;

    std::shared_ptr<AMP::LinearAlgebra::Variable> d_outputVariable;

    std::shared_ptr<AMP::LinearAlgebra::Matrix> d_matrix; // The matrix shared pointer

private:
};
} // namespace Operator
} // namespace AMP

#endif
