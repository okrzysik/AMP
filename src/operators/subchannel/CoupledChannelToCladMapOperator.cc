#include "operators/subchannel/CoupledChannelToCladMapOperator.h"
#include "operators/subchannel/CoupledChannelToCladMapOperatorParameters.h"
#include "operators/subchannel/SubchannelConstants.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {


CoupledChannelToCladMapOperator::CoupledChannelToCladMapOperator(
    const AMP::shared_ptr<CoupledChannelToCladMapOperatorParameters> &params )
    : Operator( params )
{
    AMP_ASSERT( params->d_thermalMapOperator.get() != nullptr );
    AMP_ASSERT( params->d_densityMapOperator.get() != nullptr );
    d_thermalMapOperator     = params->d_thermalMapOperator;
    d_densityMapOperator     = params->d_densityMapOperator;
    d_flowVariable           = params->d_variable;
    d_Mesh                   = params->d_subchannelMesh;
    d_subchannelPhysicsModel = params->d_subchannelPhysicsModel;
    d_subchannelTemperature  = params->d_vector;
    if ( d_Mesh.get() != nullptr ) {
        AMP_ASSERT( d_subchannelPhysicsModel != nullptr );
    }
    if ( d_subchannelTemperature.get() != nullptr ) {
        d_subchannelDensity = d_subchannelTemperature->cloneVector();
        AMP::LinearAlgebra::Variable::shared_ptr densityVariable(
            new AMP::LinearAlgebra::Variable( "Density" ) );
        d_subchannelDensity->setVariable( densityVariable );
    }
}


void CoupledChannelToCladMapOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                             AMP::LinearAlgebra::Vector::shared_ptr )
{

    // Compute the subchannel temperature and density
    if ( d_Mesh != nullptr ) {
        const double h_scale =
            1.0 /
            Subchannel::scaleEnthalpy; // Scale to change the input vector back to correct units
        const double P_scale =
            1.0 /
            Subchannel::scalePressure; // Scale to change the input vector back to correct units

        AMP::LinearAlgebra::Vector::const_shared_ptr uInternal = subsetInputVector( u );

        AMP::Discretization::DOFManager::shared_ptr faceDOFManager = uInternal->getDOFManager();
        AMP::Discretization::DOFManager::shared_ptr scalarFaceDOFManager =
            d_subchannelTemperature->getDOFManager();

        AMP::Mesh::MeshIterator face =
            AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( d_Mesh, 0 );
        // AMP::Mesh::MeshIterator end_face = face.end();

        std::vector<size_t> dofs;
        std::vector<size_t> scalarDofs;
        for ( size_t i = 0; i < face.size(); i++ ) {
            faceDOFManager->getDOFs( face->globalID(), dofs );
            scalarFaceDOFManager->getDOFs( face->globalID(), scalarDofs );
            std::map<std::string, AMP::shared_ptr<std::vector<double>>> temperatureArgMap;
            temperatureArgMap.insert(
                std::make_pair( std::string( "enthalpy" ),
                                AMP::shared_ptr<std::vector<double>>( new std::vector<double>(
                                    1, h_scale * uInternal->getValueByGlobalID( dofs[0] ) ) ) ) );
            temperatureArgMap.insert(
                std::make_pair( std::string( "pressure" ),
                                AMP::shared_ptr<std::vector<double>>( new std::vector<double>(
                                    1, P_scale * uInternal->getValueByGlobalID( dofs[1] ) ) ) ) );
            std::vector<double> temperatureResult( 1 );
            std::vector<double> specificVolume( 1 );
            d_subchannelPhysicsModel->getProperty(
                "Temperature", temperatureResult, temperatureArgMap );
            d_subchannelPhysicsModel->getProperty(
                "SpecificGeomType::Volume", specificVolume, temperatureArgMap );
            d_subchannelTemperature->setValueByGlobalID( scalarDofs[0], temperatureResult[0] );
            d_subchannelDensity->setValueByGlobalID( scalarDofs[0], 1.0 / specificVolume[0] );
            ++face;
        }

        d_subchannelTemperature->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
        d_subchannelDensity->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    }

    // Map the temperature and density back to the clad
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    d_thermalMapOperator->apply( d_subchannelTemperature, nullVec );
    d_densityMapOperator->apply( d_subchannelDensity, nullVec );
}
}
}
