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
    explicit NullOperator(const AMP::shared_ptr<OperatorParameters> & params): Operator(params) {}

    //! Destructor
    virtual ~NullOperator() { }

    //! Empty apply call
    virtual void apply( AMP::LinearAlgebra::Vector::const_shared_ptr, 
			AMP::LinearAlgebra::Vector::shared_ptr) override { }

protected :

private :

};


}
}

#endif


