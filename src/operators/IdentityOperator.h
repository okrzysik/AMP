#ifndef included_AMP_IdentityOperator
#define included_AMP_IdentityOperator

#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorParameters.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace Operator {


/**
 * Class IdentityOperator is the identity operator A(u) = u
 */
class IdentityOperator : public LinearOperator
{
public:
    //! Constructor. This resets the matrix shared pointer.
    IdentityOperator();

    /**
     * Constructor. This resets the matrix shared pointer.
     * @param [in] params
     */
    explicit IdentityOperator( const std::shared_ptr<OperatorParameters> &params );

    //! Destructor
    virtual ~IdentityOperator() {}

    /**
     * The apply function for this operator, A, performs the following operation:
     * r = A(u)
     * Here, A(u) is simply a Matrix-Vector multiplication.
     * @param [in] u input vector.
     * @param [out] f output vector.
     */
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    /**
     * This function is useful for re-initializing/updating an operator
     * \param params
     *    parameter object containing parameters to change
     */
    void reset( const std::shared_ptr<OperatorParameters> &params ) override;

    /**
     * Copies the shared pointer for the matrix representation of this linear operator.
     *  @param [in] in_mat The matrix representation of this linear operator.
     */
    void setMatrix( std::shared_ptr<AMP::LinearAlgebra::Matrix> in_mat ) override;

    //! Set the input variable
    virtual void setInputVariable( AMP::LinearAlgebra::Variable::shared_ptr var )
    {
        d_inputVariable = var;
    }

    //! Set the output variable
    virtual void setOutputVariable( AMP::LinearAlgebra::Variable::shared_ptr var )
    {
        d_outputVariable = var;
    }

private:
};
} // namespace Operator
} // namespace AMP

#endif
