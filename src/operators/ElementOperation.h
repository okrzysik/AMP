
#ifndef included_AMP_ElementOperation
#define included_AMP_ElementOperation

#include "AMP/operators/ElementOperationParameters.h"

#include <memory>


namespace AMP::Operator {

/**
  An abstract base class for representing the element level computation
  performed within a finite element operator. Concrete implementations
  must implement the apply() function.
  */
class ElementOperation
{
public:
    /**
      Constructor.
      */
    explicit ElementOperation( std::shared_ptr<const ElementOperationParameters> ) {}

    /**
      Destructor.
      */
    virtual ~ElementOperation() {}

    /**
      This is where the element level computation in a FE operator is performed. Each derived
      class must provide an implementation that is appropriate for use within its respective FE
      operator.
      */
    virtual void apply() = 0;

protected:
private:
};
} // namespace AMP::Operator

#endif
