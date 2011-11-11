#include "ElementOperationFactory.h"
#include "utils/Utilities.h"
#include "MassLinearElement.h"
#include "SourceNonlinearElement.h"
//#include "ConsMomentumGalWFLinearElement.h"
//#include "ConsMassGalWFLinearElement.h"
//#include "MechanicsLinearElement.h"
//#include "MechanicsNonlinearElement.h"
//#include "MechanicsLinearUpdatedLagrangianElement.h"
//#include "MechanicsNonlinearUpdatedLagrangianElement.h"
//#include "DiffusionLinearElement.h"
//#include "DiffusionNonlinearElement.h"
//#include "operators/mechanics/MechanicsElement.h"
//#include "operators/diffusion/DiffusionElement.h"

namespace AMP {
namespace Operator {

boost::shared_ptr<ElementOperation>
ElementOperationFactory::createElementOperation(boost::shared_ptr<Database>  elementOperationDb )
  {
AMP_ERROR("ElementOperationFactory is not converted yet");
/*
    boost::shared_ptr<ElementOperation> retElementOp;
    boost::shared_ptr<ElementOperationParameters> params;
    
    AMP_INSIST(elementOperationDb.get()!=NULL, "ElementOperationFactory::createElementOperation:: NULL Database object input");

    std::string name = elementOperationDb->getString("name");
    
    params.reset( new ElementOperationParameters(elementOperationDb));

    if(name=="MechanicsLinearElement")
      {
    retElementOp.reset(new MechanicsLinearElement(params));
      }
    else if(name=="MechanicsNonlinearElement")
      {
    retElementOp.reset(new MechanicsNonlinearElement(params));
      }
    else if(name=="MechanicsLinearUpdatedLagrangianElement")
      {
    retElementOp.reset(new MechanicsLinearUpdatedLagrangianElement(params));
      }
    else if(name=="MechanicsNonlinearUpdatedLagrangianElement")
      {
    retElementOp.reset(new MechanicsNonlinearUpdatedLagrangianElement(params));
      }
    else if (name=="DiffusionLinearElement")
      {
    retElementOp.reset(new DiffusionLinearElement(params));
      }
    else if (name=="DiffusionNonlinearElement")
      {
    retElementOp.reset(new DiffusionNonlinearElement(params));
      }
    else if (name=="ConsMomentumGalWFLinearElement")
      {
    retElementOp.reset(new ConsMomentumGalWFLinearElement(params));
      }
    else if (name=="ConsMassGalWFLinearElement")
      {
    retElementOp.reset(new ConsMassGalWFLinearElement(params));
      }
    else if (name=="MassLinearElement")
      {
    retElementOp.reset(new MassLinearElement(params));
      }
    else if (name=="SourceNonlinearElement")
      {
    retElementOp.reset(new SourceNonlinearElement(params));
      }
    
    return retElementOp;
*/
}
}
  
}
