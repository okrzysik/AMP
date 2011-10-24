
#ifndef included_AMP_ColumnOperatorParameters
#define included_AMP_ColumnOperatorParameters

/*AMP files */
#include "operators/OperatorParameters.h"

/*Boost files */
#include "boost/shared_ptr.hpp"

#include <vector>

namespace AMP {
namespace Operator {

  /**
    A class that encapsulates the parameters required to construct
    the composite Operator operator.
    @see ColumnOperator
    */
  class ColumnOperatorParameters : public OperatorParameters {
    public :

      ColumnOperatorParameters(const boost::shared_ptr<AMP::Database>& db)
        : OperatorParameters(db) { }

      virtual ~ColumnOperatorParameters() { }

      std::vector< boost::shared_ptr< OperatorParameters > > d_OperatorParameters;
  };

}  
}


#endif



