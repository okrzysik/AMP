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
    AMP_ASSERT(params->d_thermalMapOperator.get()!=NULL);
    AMP_ASSERT(params->d_densityMapOperator.get()!=NULL);
    d_thermalMapOperator            = params->d_thermalMapOperator; 
    d_densityMapOperator            = params->d_densityMapOperator; 
    d_flowVariable           = params->d_variable; 
    d_Mesh                   = params->d_subchannelMesh; 
    d_subchannelPhysicsModel = params->d_subchannelPhysicsModel; 
    d_subchannelTemperature = params->d_vector; 
    if ( d_subchannelTemperature.get()!=NULL ) {
        d_subchannelDensity = d_subchannelTemperature->cloneVector();
        AMP::LinearAlgebra::Variable::shared_ptr densityVariable( new AMP::LinearAlgebra::Variable("Density") );
        d_subchannelDensity->setVariable( densityVariable );
    }
}


void CoupledChannelToCladMapOperator :: apply( AMP::LinearAlgebra::Vector::const_shared_ptr ,
    AMP::LinearAlgebra::Vector::const_shared_ptr u, AMP::LinearAlgebra::Vector::shared_ptr ,
    const double a, const double b)
{

    // Compute the subchannel temperature and density
    if ( d_Mesh!=NULL) {
        const double h_scale = 1.0/Subchannel::scaleEnthalpy;                 // Scale to change the input vector back to correct units
        const double P_scale = 1.0/Subchannel::scaleEnthalpy;                 // Scale to change the input vector back to correct units

        AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );

        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = uInternal->getDOFManager(); 
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager = d_subchannelTemperature->getDOFManager(); 

        AMP::Mesh::MeshIterator face = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(d_Mesh, 0);
        AMP::Mesh::MeshIterator end_face = face.end();
        
        std::vector<size_t> dofs;
        std::vector<size_t> scalarDofs;
        for(size_t i=0; i<face.size(); i++){
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, boost::shared_ptr<std::vector<double> > > temperatureArgMap;
            temperatureArgMap.insert(std::make_pair("enthalpy",new std::vector<double>(1,h_scale*uInternal->getValueByGlobalID(dofs[0]))));
            temperatureArgMap.insert(std::make_pair("pressure",new std::vector<double>(1,P_scale*uInternal->getValueByGlobalID(dofs[1]))));
            std::vector<double> temperatureResult(1);
            std::vector<double> specificVolume(1);
            d_subchannelPhysicsModel->getProperty("Temperature", temperatureResult, temperatureArgMap); 
            d_subchannelPhysicsModel->getProperty("SpecificVolume", specificVolume, temperatureArgMap);
            d_subchannelTemperature->setValueByGlobalID(scalarDofs[0], temperatureResult[0]);
            d_subchannelDensity->setValueByGlobalID(scalarDofs[0],1.0/specificVolume[0]);
            ++face;
        }
        
        d_subchannelTemperature->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        d_subchannelDensity->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        
    }

    // Map the temperature and density back to the clad
    AMP::LinearAlgebra::Vector::shared_ptr   nullVec;
    d_thermalMapOperator->apply(nullVec, d_subchannelTemperature, nullVec, a, b);
    d_densityMapOperator->apply(nullVec, d_subchannelDensity, nullVec, a, b);

}


}
}



