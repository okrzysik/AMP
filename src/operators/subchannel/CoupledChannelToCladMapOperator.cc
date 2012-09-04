#include "CoupledChannelToCladMapOperator.h"
#include "CoupledChannelToCladMapOperatorParameters.h"
#include "operators/subchannel/SubchannelConstants.h"

#include "utils/Utilities.h"
#include "ampmesh/StructuredMeshHelper.h"


namespace AMP {
namespace Operator {


CoupledChannelToCladMapOperator::CoupledChannelToCladMapOperator(const boost::shared_ptr<CoupledChannelToCladMapOperatorParameters>& params)
      : Operator(params)
{
    d_mapOperator            = params->d_mapOperator; 
    d_flowVariable           = params->d_variable; 
    d_Mesh                   = params->d_subchannelMesh; 
    d_subchannelPhysicsModel = params->d_subchannelPhysicsModel; 
    d_subchannelTemperature = params->d_vector; 
    // AMP::LinearAlgebra::Vector::shared_ptr temperature  = params->d_vector ; 
    // AMP::LinearAlgebra::VS_Mesh meshSelector("meshSelector", d_Mesh);
    // d_subchannelTemperature = temperature->select(meshSelector, d_mapOperator->getOutputVariable()->getName());
}


void CoupledChannelToCladMapOperator :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr ,
    AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr ,
    const double a, const double b)
{

    // Compute the subchannel temperature       
    if ( d_Mesh!=NULL) {
        const double h_scale = 1.0/Subchannel::scaleEnthalpy;                 // Scale to change the input vector back to correct units
        const double P_scale = 1.0/Subchannel::scaleEnthalpy;                 // Scale to change the input vector back to correct units

        AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );

        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = uInternal->getDOFManager(); 
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = d_subchannelTemperature->getDOFManager(); 

        AMP::Mesh::MeshIterator face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(d_Mesh, 0);
        AMP::Mesh::MeshIterator end_face = face.end();
        
        for( ; face != end_face; ++face){
            std::vector<size_t> dofs;
            std::vector<size_t> scalarDofs;
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
            temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_scale*uInternal->getValueByGlobalID(dofs[0]))));
            temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,P_scale*uInternal->getValueByGlobalID(dofs[1]))));
            std::vector<double> temperatureResult(1);
            d_subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
            d_subchannelTemperature->setValueByGlobalID(scalarDofs[0], temperatureResult[0]);
        }
        
        d_subchannelTemperature->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        
    }

    // Map the temperature back to the clad
    AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
    d_mapOperator->apply(nullVec, d_subchannelTemperature, nullVec, a, b);

}


}
}



