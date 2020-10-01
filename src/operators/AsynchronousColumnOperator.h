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
    /** Constructor
     */
    explicit AsynchronousColumnOperator( const std::shared_ptr<OperatorParameters> & );

    void apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
		AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f );

    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr f );

    void append( std::shared_ptr<Operator> op ) override;
};
} // namespace Operator
} // namespace AMP

#endif
