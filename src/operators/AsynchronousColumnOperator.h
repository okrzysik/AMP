#ifndef included_AMP_AsynchronousColumnOperator
#define included_AMP_AsynchronousColumnOperator

#include "ColumnOperator.h"

namespace AMP {
namespace Operator {

/** \brief  A column operator of asynchronous operators.  The apply method will start the list
    of operators then finalize the list of operators
    */
class AsynchronousColumnOperator : public ColumnOperator
{
public:
    //! Empty constructor
    explicit AsynchronousColumnOperator();

    //! Default constructor
    explicit AsynchronousColumnOperator( std::shared_ptr<const OperatorParameters> );

    //! Return the name of the operator
    std::string type() const override { return "AsynchronousColumnOperator"; }

    //! The apply routine for the column operator calls apply on each of the component operators
    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    //! Start an apply operation
    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr f );

    //! Finish an apply operation (arguments should match applyStart)
    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f );

    // Append the given operator
    void append( std::shared_ptr<Operator> op ) override;
};
} // namespace Operator
} // namespace AMP

#endif
