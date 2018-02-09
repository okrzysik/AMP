
#ifndef included_AMP_ElementOperation
#define included_AMP_ElementOperation

#include "ElementOperationParameters.h"

#include "AMP/utils/shared_ptr.h"

namespace AMP {
namespace Operator {

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
    explicit ElementOperation( const AMP::shared_ptr<ElementOperationParameters> & ) {}

    /**
      Destructor.
      */
    virtual ~ElementOperation() {}

    /**
      This is where the element level computation in a FE operator is performed. Each derived
      class must provide an implementation that is appropriate for use within its respective FE
      operator.
      */
    virtual void apply()
    {
        // Implemented in derived classes.
    }

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
