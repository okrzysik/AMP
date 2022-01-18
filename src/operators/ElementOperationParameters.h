
#ifndef included_AMP_ElementOperationParameters
#define included_AMP_ElementOperationParameters

#include <memory>

#include "AMP/utils/Database.h"

#include "AMP/utils/ParameterBase.h"

namespace AMP::Operator {

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
    explicit ElementOperationParameters( std::shared_ptr<AMP::Database> db ) : d_db( db ) {}

    /**
      Destructor.
      */
    virtual ~ElementOperationParameters() {}

    std::shared_ptr<AMP::Database> d_db; /**< Database object which needs to be
                                             initialized specific to the element operation. */

protected:
private:
};
} // namespace AMP::Operator

#endif
