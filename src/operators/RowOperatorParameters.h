#ifndef included_RowOperatorParameters
#define included_RowOperatorParameters

#include "utils/shared_ptr.h"
#include "operators/OperatorParameters.h"
#include "operators/Operator.h"

namespace AMP{
namespace Operator {

class RowOperatorParameters: public OperatorParameters
{
  public:
  
    RowOperatorParameters(const AMP::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) { }

    virtual ~RowOperatorParameters(){ };


    std::vector< AMP::shared_ptr<AMP::Operator> > d_Operator;

    std::vector< AMP::shared_ptr<AMP::OperatorParameters> > d_OperatorParameters;

    std::vector< double > scalea;

};


}  
}

#endif

