
#ifndef included_AMP_BlockOperatorParameters
#define included_AMP_BlockOperatorParameters

#include "operators/OperatorParameters.h"

#include "boost/shared_ptr.hpp"

namespace AMP {
  namespace Operator {

    class BlockOperatorParameters : public OperatorParameters {
      public :

        BlockOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
          : OperatorParameters(db) { }

        virtual ~BlockOperatorParameters() { }

        std::vector<std::vector<boost::shared_ptr<OperatorParameters> > > d_blockParams;
    };

  }  
}

#endif



