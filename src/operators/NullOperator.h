#ifndef included_AMP_NullOperator
#define included_AMP_NullOperator

#include "utils/shared_ptr.h"

#include "operators/Operator.h"
#include "operators/OperatorParameters.h"

#include "vectors/Vector.h"

#include <string>



namespace AMP {
namespace Operator {


/**
  * Class NullOperator is an empty operator that does nothing
  */
class NullOperator : public AMP::Operator::Operator {
public :

    typedef AMP::shared_ptr<AMP::Operator::Operator>  shared_ptr;

    //! Default constructor
    NullOperator(void) {}

    //! Constructor
    NullOperator(const AMP::shared_ptr<OperatorParameters> & params): Operator(params) {}

    //! Destructor
    virtual ~NullOperator() { }

    //! Empty apply call
    virtual void apply(AMP::LinearAlgebra::Vector::const_shared_ptr f, 
        AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr r,
        const double a = -1.0, const double b = 1.0) {}

protected :

private :

};


}
}

#endif


