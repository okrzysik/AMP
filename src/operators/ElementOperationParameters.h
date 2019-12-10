
#ifndef included_AMP_ElementOperationParameters
#define included_AMP_ElementOperationParameters

#include "AMP/utils/shared_ptr.h"

#include "AMP/utils/Database.h"

#include "AMP/utils/ParameterBase.h"

namespace AMP {
namespace Operator {

/**
 An abstract base class that encapsulates parameters used to initialize the ElementOperation
 used within a FEOperator.
 @see ElementOperation
 */
class ElementOperationParameters : public ParameterBase
{
public:
    /**
      Constructor.
      */
    explicit ElementOperationParameters( AMP::shared_ptr<AMP::Database> db ) : d_db( db ) {}

    /**
      Destructor.
      */
    virtual ~ElementOperationParameters() {}

    AMP::shared_ptr<AMP::Database> d_db; /**< Database object which needs to be
                                             initialized specific to the element operation. */

protected:
private:
};
} // namespace Operator
} // namespace AMP

#endif
