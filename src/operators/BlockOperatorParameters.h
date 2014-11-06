
#ifndef included_AMP_BlockOperatorParameters
#define included_AMP_BlockOperatorParameters

#include "operators/OperatorParameters.h"

#include "utils/shared_ptr.h"

namespace AMP {
  namespace Operator {

    class BlockOperatorParameters : public OperatorParameters {
      public :

        BlockOperatorParameters(const AMP::shared_ptr<AMP::Database> &db)
          : OperatorParameters(db) { }

        virtual ~BlockOperatorParameters() { }

        std::vector<std::vector<AMP::shared_ptr<OperatorParameters> > > d_blockParams;
    };

  }  
}

#endif



