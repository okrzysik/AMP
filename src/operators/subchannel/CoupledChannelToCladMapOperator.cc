#include "CoupledChannelToCladMapOperator.h"
#include "CoupledChannelToCladMapOperatorParameters.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {

    CoupledChannelToCladMapOperator::CoupledChannelToCladMapOperator(const boost::shared_ptr<OperatorParameters>& params)
      : Operator(params)
    {

      boost::shared_ptr<CoupledChannelToCladMapOperatorParameters> myparams = boost::dynamic_pointer_cast<CoupledChannelToCladMapOperatorParameters>(params);
      // do nothing
      d_mapOperator            = myparams->d_mapOperator ; 
      d_flowVariable           = myparams->d_variable ; 
      d_Mesh                   = myparams->d_subchannelMesh ; 
      d_subchannelTemperature  = myparams->d_vector ; 
      d_subchannelPhysicsModel = myparams->d_subchannelPhysicsModel; 

      AMP::LinearAlgebra::VS_Mesh meshSelector("meshSelector", d_Mesh);
      d_subchannelTemperature = u->select(meshSelector, d_MapOperator->getOutputVariable()->getName());

    }

    void
      CoupledChannelToCladMapOperator :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr ,
          AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr ,
          const double a, const double b)
      {
        AMP::LinearAlgebra::Vector::shared_ptr   nullVec;

        AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );

        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = uInternal->getDofManager(); 
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = d_subchannelTemperature->getDofManager(); 

        AMP::Mesh::MeshIterator face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(d_Mesh, 0);
        AMP::Mesh::MeshIterator end_face   = face.end();
        
        for( ; face != end_face; ++face){
          std::vector<size_t> dofs;
          std::vector<size_t> scalarDofs;
          faceDOFManager->getDOFs( face->globalID(), dofs );
          scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
          std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
          temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,uInternal->getValueByGlobalID(dofs[0]))));
          temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,uInternal->getValueByGlobalID(dofs[1]))));
          std::vector<double> temperatureResult(1);
          d_subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
          double compTin = temperatureResult[0];
          d_subchannelTemperature->setValueByGlobalID(dofs[1], temperatureResult[0])
        }
        
        d_subchannelTemperature->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        d_MapOperator->apply(nullVec, d_subchannelTemperature, nullVec, a, b);

      }


}
}



