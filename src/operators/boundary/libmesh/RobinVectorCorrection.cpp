#include "RobinVectorCorrection.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "ProfilerApp.h"
#include "RobinMatrixCorrectionParameters.h"

// Libmesh headers
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

#include <string>

namespace AMP {
namespace Operator {


RobinVectorCorrection::RobinVectorCorrection(
    std::shared_ptr<const NeumannVectorCorrectionParameters> params )
    : NeumannVectorCorrection( params )
{
    d_hef        = 0;
    d_alpha      = 0;
    d_beta       = 0;
    d_gamma      = 0;
    d_skipParams = false;
    reset( params );
    d_InstanceID = d_iInstance_id;
}


void RobinVectorCorrection::reset( std::shared_ptr<const OperatorParameters> params )
{
    NeumannVectorCorrection::reset( params );

    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameters" );
    AMP_INSIST( ( ( ( params->d_db ).get() ) != nullptr ), "NULL database" );

    d_skipParams = params->d_db->getWithDefault( "skip_params", false );

    AMP_INSIST( params->d_db->keyExists( "alpha" ), "Missing key: prefactor alpha" );
    d_alpha = params->d_db->getScalar<double>( "alpha" );

    AMP_INSIST( params->d_db->keyExists( "beta" ), "Missing key: prefactor beta" );
    d_beta = params->d_db->getScalar<double>( "beta" );

    AMP_INSIST( params->d_db->keyExists( "gamma" ), "Missing key: total prefactor gamma" );
    d_gamma = params->d_db->getScalar<double>( "gamma" );
}


void RobinVectorCorrection::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                   AMP::LinearAlgebra::Vector::shared_ptr r )
{
    PROFILE_START( "apply" );
    AMP_INSIST( ( ( r.get() ) != nullptr ), "NULL Residual Vector" );
    AMP_INSIST( ( ( u.get() ) != nullptr ), "NULL Solution Vector" );

    auto rInternal = this->subsetInputVector( r );
    auto uInternal = this->subsetInputVector( u );

    AMP_ASSERT( uInternal->getUpdateStatus() ==
                AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
    // rInternal->makeConsistent ( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    std::vector<std::string> variableNames;
    size_t numVar = 0;
    if ( d_robinPhysicsModel ) {
        variableNames = d_robinPhysicsModel->getVariableName();
        numVar        = variableNames.size();
    }

    d_elementInputVec.resize( numVar + 1 );
    d_elementInputVec[0] = d_variableFlux;

    if ( d_robinPhysicsModel ) {
        for ( size_t i = 0; i < variableNames.size(); i++ ) {
            std::string cview = variableNames[i] + " view";
            if ( d_Frozen ) {
                if ( d_Frozen->select( AMP::LinearAlgebra::VS_ByVariableName( variableNames[i] ),
                                       cview ) != nullptr ) {
                    d_elementInputVec[i + 1] = d_Frozen->constSelect(
                        AMP::LinearAlgebra::VS_ByVariableName( variableNames[i] ), cview );
                } else {
                    d_elementInputVec[i + 1] = uInternal->constSelect(
                        AMP::LinearAlgebra::VS_ByVariableName( variableNames[i] ), cview );
                }
            } else {
                d_elementInputVec[i + 1] = uInternal->constSelect(
                    AMP::LinearAlgebra::VS_ByVariableName( variableNames[i] ), cview );
            }
            AMP_INSIST( d_elementInputVec[i + 1],
                        "Did not find vector '" + variableNames[i] + "'" );
            AMP_ASSERT( d_elementInputVec[i + 1]->getUpdateStatus() ==
                        AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED );
        }

        if ( d_iDebugPrintInfoLevel == 100 ) {
            std::cout << "processing robin boundary operator " << d_InstanceID << "\n";
        }
    }

    // Get the DOF managers
    auto dofManager = rInternal->getDOFManager();
    std::shared_ptr<AMP::Discretization::DOFManager> gpDOFManager;
    if ( d_isFluxGaussPtVector && d_variableFlux != nullptr )
        gpDOFManager = d_variableFlux->getDOFManager();

    // Check that the DOF managers match for the different vectors
    AMP_ASSERT( *dofManager == *( uInternal->getDOFManager() ) );
    AMP_ASSERT( *dofManager == *( rInternal->getDOFManager() ) );
    if ( !d_isFluxGaussPtVector ) {
        if ( d_variableFlux )
            AMP_ASSERT( *dofManager == *( d_variableFlux->getDOFManager() ) );
    }

    auto numIds = d_boundaryIds.size();
    std::vector<size_t> gpDofs;
    std::vector<size_t> dofs;
    std::vector<std::vector<size_t>> dofIndices;
    std::vector<size_t> dofsElementVec;
    PROFILE_START( "integration loop" );
    for ( unsigned int nid = 0; nid < numIds; nid++ ) {
        unsigned int numDofIds = d_dofIds[nid].size();

        for ( unsigned int k = 0; k < numDofIds; k++ ) {

            auto bnd1 =
                d_Mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, d_boundaryIds[nid], 0 );

            for ( const auto &elem : bnd1 ) {
                PROFILE_START( "prepare element", 2 );

                // Get the nodes for the current element
                d_currNodes             = elem.getElements( AMP::Mesh::GeomType::Vertex );
                auto numNodesInCurrElem = d_currNodes.size();

                dofIndices.resize( numNodesInCurrElem );
                // Get the dofs for the vectors
                std::vector<AMP::Mesh::MeshElementID> ids( d_currNodes.size() );
                for ( size_t i = 0; i < d_currNodes.size(); i++ )
                    ids[i] = d_currNodes[i].globalID();

                for ( unsigned int i = 0; i < numNodesInCurrElem; i++ )
                    dofManager->getDOFs( d_currNodes[i].globalID(), dofIndices[i] );

                dofs.resize( numNodesInCurrElem );
                for ( size_t n = 0; n < dofIndices.size(); n++ )
                    dofs[n] = dofIndices[n][d_dofIds[nid][k]];

                AMP_ASSERT( dofs.size() == numNodesInCurrElem );

                if ( d_isFluxGaussPtVector && d_IsCoupledBoundary[nid] ) {
                    gpDOFManager->getDOFs( elem.globalID(), gpDofs );
                    AMP_ASSERT( gpDofs.size() > 0 );
                }

                // Get the current libmesh element
                auto fe   = d_libmeshElements.getFEBase( elem.globalID() );
                auto rule = d_libmeshElements.getQBase( elem.globalID() );
                AMP_ASSERT( fe != nullptr );
                AMP_ASSERT( rule != nullptr );
                auto numGaussPts = rule->n_points();

                auto JxW = fe->get_JxW();
                auto phi = fe->get_phi();
                PROFILE_STOP( "prepare element", 2 );

                std::vector<std::vector<double>> inputArgs(
                    d_elementInputVec.size(), std::vector<double>( numNodesInCurrElem ) );
                std::vector<std::vector<double>> inputArgsAtGpts(
                    d_elementInputVec.size(), std::vector<double>( numGaussPts ) );
                std::vector<double> beta( numGaussPts, d_beta );
                std::vector<double> gamma( numGaussPts, d_gamma );
                PROFILE_START( "get conductance", 2 );
                if ( d_robinPhysicsModel ) {
                    unsigned int startIdx = 0;
                    if ( d_isFluxGaussPtVector && d_IsCoupledBoundary[nid] ) {
                        d_variableFlux->getValuesByGlobalID(
                            gpDofs.size(), &gpDofs[0], &inputArgsAtGpts[0][0] );
                        startIdx = 1;
                    }
                    for ( unsigned int m = startIdx; m < d_elementInputVec.size(); m++ ) {
                        // Note: elementInputVecs may use different DOFManagers from u and r
                        // internal
                        d_elementInputVec[m]->getDOFManager()->getDOFs( ids, dofsElementVec );
                        AMP_ASSERT( dofsElementVec.size() == dofs.size() );
                        d_elementInputVec[m]->getValuesByGlobalID(
                            dofs.size(), &dofsElementVec[0], &inputArgs[m][0] );
                        for ( size_t qp = 0; qp < numGaussPts; qp++ ) {
                            for ( size_t n = 0; n < numNodesInCurrElem; n++ ) {
                                inputArgsAtGpts[m][qp] += phi[n][qp] * inputArgs[m][n];
                            }
                        }
                    }

                    d_robinPhysicsModel->getConductance( beta, gamma, inputArgsAtGpts );
                }
                PROFILE_STOP( "get conductance", 2 );
                PROFILE_START( "perform integration", 2 );
                std::vector<double> values( dofs.size(), 0.0 );
                std::vector<double> gpValues( gpDofs.size(), 0.0 );
                std::vector<double> addValues( dofs.size(), 0.0 );
                uInternal->getValuesByGlobalID( dofs.size(), &dofs[0], &values[0] );
                for ( unsigned int qp = 0; qp < numGaussPts; qp++ ) {
                    libMesh::Real phi_val = 0.0;
                    for ( unsigned int l = 0; l < numNodesInCurrElem; l++ )
                        phi_val += phi[l][qp] * values[l];
                    for ( unsigned int j = 0; j < numNodesInCurrElem; j++ )
                        addValues[j] += ( JxW[qp] * phi[j][qp] * beta[qp] * phi_val );
                } // end for qp
                if ( d_IsCoupledBoundary[nid] ) {
                    if ( !d_isFluxGaussPtVector ) {
                        d_variableFlux->getValuesByGlobalID( dofs.size(), &dofs[0], &values[0] );
                    } else {
                        d_variableFlux->getValuesByGlobalID(
                            gpDofs.size(), &gpDofs[0], &gpValues[0] );
                    }
                    for ( unsigned int qp = 0; qp < numGaussPts; qp++ ) {
                        libMesh::Real phi_val = 0.0;
                        if ( !d_isFluxGaussPtVector ) {
                            for ( unsigned int l = 0; l < numNodesInCurrElem; l++ )
                                phi_val += phi[l][qp] * values[l];
                        } else {
                            phi_val = gpValues[qp];
                        }
                        for ( unsigned int j = 0; j < numNodesInCurrElem; j++ )
                            addValues[j] += -1 * ( JxW[qp] * phi[j][qp] * gamma[qp] * phi_val );
                    } // end for qp
                }     // coupled
                rInternal->addValuesByGlobalID( dofs.size(), &dofs[0], &addValues[0] );
                PROFILE_STOP( "perform integration", 2 );

            } // end for bnd
        }     // end for dof ids
    }         // end for nid
    PROFILE_STOP( "integration loop" );

    rInternal->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_ADD );
    // std::cout << rInternal << std::endl;

    PROFILE_STOP( "apply" );
}


std::shared_ptr<OperatorParameters>
    RobinVectorCorrection::getJacobianParameters( AMP::LinearAlgebra::Vector::const_shared_ptr )
{
    auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
    tmp_db->putScalar( "skip_params", true );
    tmp_db->putScalar( "skip_rhs_correction", true );
    tmp_db->putScalar( "skip_matrix_correction", false );
    tmp_db->putScalar( "IsFluxGaussPtVector", d_isFluxGaussPtVector );
    auto outParams                 = std::make_shared<RobinMatrixCorrectionParameters>( tmp_db );
    outParams->d_robinPhysicsModel = d_robinPhysicsModel;
    outParams->d_elementInputVec   = d_elementInputVec;
    outParams->d_variableFlux      = d_variableFlux;

    return outParams;
}
} // namespace Operator
} // namespace AMP
