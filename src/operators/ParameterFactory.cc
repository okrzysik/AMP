#include "ParameterFactory.h"
#include "DirichletMatrixCorrectionParameters.h"
#include "MechanicsLinearFEOperatorParameters.h"
#include "MechanicsNonlinearFEOperatorParameters.h"
#include "NeutronicsRhsParameters.h"

namespace AMP {
namespace Operator {

boost::shared_ptr<OperatorParameters>
ParameterFactory::createParameter(boost::shared_ptr<AMP::Database>  input_db, AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh)
{
  boost::shared_ptr<OperatorParameters> retParameters;
  std::string name;

  AMP_INSIST(input_db.get()!=NULL, "ParameterFactory::createParameter:: NULL Database object input");
  AMP_INSIST(input_db->keyExists("name"),  "ParameterFactory::createParameter:: key 'name' must be a part of database ");

  name = input_db->getString("name");
  
  if(name=="DirichletMatrixCorrection")
    {
      retParameters.reset(new DirichletMatrixCorrectionParameters(input_db));
    }
  else if (name=="MechanicsLinearFEOperator")
    {
      retParameters.reset(new MechanicsLinearFEOperatorParameters(input_db));
    }
  else if (name=="MechanicsNonlinearFEOperator")
    {
      retParameters.reset(new MechanicsNonlinearFEOperatorParameters(input_db));
    }
  else if (name=="NeutronicsRhs")
    {
      retParameters.reset(new NeutronicsRhsParameters(input_db));
    }
  else
    {
    }

  retParameters->d_MeshAdapter = mesh;
  
  return retParameters;
}
  
}
}

