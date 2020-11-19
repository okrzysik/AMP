#ifndef included_AMP_NullOperator
#define included_AMP_NullOperator

#include <memory>

#include "AMP/operators/Operator.h"
#include "AMP/operators/OperatorParameters.h"

#include "AMP/vectors/Vector.h"

#include <string>


namespace AMP {
namespace Operator {


/**
 * Class NullOperator is an empty operator that does nothing
 */
class NullOperator : public AMP::Operator::Operator
{
public:
    typedef std::shared_ptr<AMP::Operator::Operator> shared_ptr;

    //! Default constructor
    NullOperator( void ) {}

    //! Constructor
    explicit NullOperator( const std::shared_ptr<OperatorParameters> &params ) : Operator( params )
    {
    }

    //! Destructor
    virtual ~NullOperator() {}

    //! Return the name of the operator
    std::string type() const override { return "NullOperator"; }

    //! Empty apply call
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr,
                        AMP::LinearAlgebra::Vector::shared_ptr ) override
    {
    }

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
