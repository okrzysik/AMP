#ifndef included_RowOperatorParameters
#define included_RowOperatorParameters

#include "boost/shared_ptr.hpp"
#include "operators/OperatorParameters.h"
#include "operators/Operator.h"

namespace AMP{
namespace Operator {

class RowOperatorParameters: public OperatorParameters
{
  public:
  
    RowOperatorParameters(const boost::shared_ptr<AMP::Database> &db)
        : OperatorParameters(db) { }

    virtual ~RowOperatorParameters(){ };


    std::vector< boost::shared_ptr<AMP::Operator> > d_Operator;

    std::vector< boost::shared_ptr<AMP::OperatorParameters> > d_OperatorParameters;

    std::vector< double > scalea;

};


}  
}

#endif

