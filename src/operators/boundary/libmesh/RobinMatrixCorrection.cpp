#include "RobinMatrixCorrection.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "RobinMatrixCorrectionParameters.h"

// Libmesh files
DISABLE_WARNINGS
#include "libmesh/auto_ptr.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/face_quad4.h"
#include "libmesh/node.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

#include <string>

namespace AMP::Operator {

RobinMatrixCorrection::RobinMatrixCorrection(
    std::shared_ptr<const RobinMatrixCorrectionParameters> params )
    : BoundaryOperator( params )
{
    d_hef   = 0;
    d_alpha = 0;
    d_beta  = 0;
    d_gamma = 0;
    d_JxW   = nullptr;
    d_phi   = nullptr;

    auto feTypeOrderName = params->d_db->getWithDefault<std::string>( "FE_ORDER", "FIRST" );
    auto feFamilyName    = params->d_db->getWithDefault<std::string>( "FE_FAMILY", "LAGRANGE" );
    auto qruleTypeName   = params->d_db->getWithDefault<std::string>( "QRULE_TYPE", "QGAUSS" );
    d_qruleOrderName     = params->d_db->getWithDefault<std::string>( "QRULE_ORDER", "DEFAULT" );

    d_feTypeOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>( feTypeOrderName );
    d_feFamily    = libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( feFamilyName );
    d_qruleType   = libMesh::Utility::string_to_enum<libMeshEnums::QuadratureType>( qruleTypeName );

    d_variable = params->d_variable;

    d_NeumannParams.reset( new AMP::Operator::NeumannVectorCorrectionParameters( params->d_db ) );
    d_NeumannParams->d_variable     = params->d_variable;
    d_NeumannParams->d_Mesh         = params->d_Mesh;
    d_NeumannParams->d_variableFlux = params->d_variableFlux;
    d_NeumannCorrection.reset( new NeumannVectorCorrection( d_NeumannParams ) );

    reset( params );
}

void RobinMatrixCorrection::reset( std::shared_ptr<const OperatorParameters> params )
{

    auto myparams = std::dynamic_pointer_cast<const RobinMatrixCorrectionParameters>( params );

    AMP_INSIST( myparams.get(), "NULL parameters" );
    AMP_INSIST( myparams->d_db, "NULL database" );

    AMP_INSIST( myparams->d_db, "NULL database" );
    bool skipParams = myparams->d_db->getWithDefault<bool>( "skip_params", true );

    bool d_isFluxGaussPtVector =
        myparams->d_db->getWithDefault<bool>( "IsFluxGaussPtVector", true );

    if ( !skipParams ) {
        AMP_INSIST( myparams->d_db->keyExists( "alpha" ), "Missing key: prefactor alpha" );
        d_alpha = myparams->d_db->getScalar<double>( "alpha" );
        AMP_INSIST( d_alpha != 0.0, "prefactor alpha must be != 0.0" );

        AMP_INSIST( myparams->d_db->keyExists( "beta" ), "Missing key: prefactor beta" );
        d_beta = myparams->d_db->getScalar<double>( "beta" );

        AMP_INSIST( myparams->d_db->keyExists( "gamma" ), "Missing key: total prefactor gamma" );
        d_gamma = myparams->d_db->getScalar<double>( "gamma" );

        AMP_INSIST( params->d_db->keyExists( "number_of_ids" ),
                    "Key ''number_of_ids'' is missing!" );
        int numIds = params->d_db->getScalar<int>( "number_of_ids" );

        d_boundaryIds.resize( numIds );
        d_dofIds.resize( numIds );
        d_robinValues.resize( numIds );

        char key[100];
        for ( int j = 0; j < numIds; j++ ) {

            snprintf( key, sizeof key, "id_%d", j );
            AMP_INSIST( myparams->d_db->keyExists( key ), "Key is missing!" );
            d_boundaryIds[j] = myparams->d_db->getScalar<int>( key );

            snprintf( key, sizeof key, "number_of_dofs_%d", j );
            AMP_INSIST( myparams->d_db->keyExists( key ), "Key is missing!" );
            int numDofIds = myparams->d_db->getScalar<int>( key );

            d_dofIds[j].resize( numDofIds );
            d_robinValues[j].resize( numDofIds );
            for ( int i = 0; i < numDofIds; i++ ) {
                snprintf( key, sizeof key, "dof_%d_%d", j, i );
                AMP_INSIST( myparams->d_db->keyExists( key ), "Key is missing!" );
                d_dofIds[j][i] = myparams->d_db->getScalar<int>( key );

                snprintf( key, sizeof key, "value_%d_%d", j, i );
                AMP_INSIST( myparams->d_db->keyExists( key ), "Key is missing!" );
                d_robinValues[j][i] = myparams->d_db->getScalar<double>( key );
            }
        }
    }

    d_robinPhysicsModel = myparams->d_robinPhysicsModel;

    d_NeumannParams->d_db->putScalar( "constant_flux",
                                      myparams->d_db->getWithDefault<bool>( "constant_flux", true ),
                                      {},
                                      AMP::Database::Check::Overwrite );
    d_NeumannParams->d_variableFlux      = myparams->d_variableFlux;
    d_NeumannParams->d_robinPhysicsModel = myparams->d_robinPhysicsModel;
    d_NeumannParams->d_db->putScalar( "gamma", d_gamma, {}, AMP::Database::Check::Overwrite );
    d_NeumannCorrection->reset( d_NeumannParams );

    bool skipMatrixCorrection =
        myparams->d_db->getWithDefault<bool>( "skip_matrix_correction", false );
    if ( !skipMatrixCorrection ) {
        // Create the libmesh elements
        AMP::Mesh::MeshIterator iterator;
        for ( auto &elem : d_boundaryIds ) {
            auto iterator2 = d_Mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, elem, 0 );
            iterator = AMP::Mesh::Mesh::getIterator( AMP::Mesh::SetOP::Union, iterator, iterator2 );
        }
        libmeshElements.reinit( iterator );

        auto inputMatrix = myparams->d_inputMatrix;
        AMP_INSIST( ( ( inputMatrix.get() ) != nullptr ), "NULL matrix" );

        auto inVec   = inputMatrix->getRightVector();
        d_dofManager = inVec->getDOFManager();

        unsigned int numIds  = d_boundaryIds.size();
        auto elementInputVec = myparams->d_elementInputVec;

        std::vector<size_t> gpDofs, dofsElementVec;
        std::vector<std::vector<size_t>> dofIndices;
        std::vector<size_t> dofs;

        std::shared_ptr<AMP::Discretization::DOFManager> gpDOFManager;
        if ( d_isFluxGaussPtVector && myparams->d_variableFlux )
            gpDOFManager = myparams->d_variableFlux->getDOFManager();

        for ( unsigned int nid = 0; nid < numIds; nid++ ) {
            unsigned int numDofIds = d_dofIds[nid].size();

            for ( unsigned int k = 0; k < numDofIds; k++ ) {
                auto bnd1 = d_Mesh->getBoundaryIDIterator(
                    AMP::Mesh::GeomType::Face, d_boundaryIds[nid], 0 );
                auto end_bnd1 = bnd1.end();
                for ( ; bnd1 != end_bnd1; ++bnd1 ) {

                    auto d_feType = std::make_shared<libMesh::FEType>( d_feTypeOrder, d_feFamily );
                    std::shared_ptr<libMesh::FEBase> d_fe(
                        ( libMesh::FEBase::build( 2, ( *d_feType ) ) ).release() );

                    if ( d_qruleOrderName == "DEFAULT" ) {
                        d_qruleOrder = d_feType->default_quadrature_order();
                    } else {
                        d_qruleOrder = libMesh::Utility::string_to_enum<libMeshEnums::Order>(
                            d_qruleOrderName );
                    }
                    std::shared_ptr<libMesh::QBase> d_qrule(
                        ( libMesh::QBase::build( d_qruleType, 2, d_qruleOrder ) ).release() );

                    // Get the nodes for the element and their global ids
                    auto currNodes = bnd1->getElements( AMP::Mesh::GeomType::Vertex );
                    dofIndices.resize( currNodes.size() );
                    std::vector<AMP::Mesh::MeshElementID> globalIDs( currNodes.size() );
                    for ( size_t j = 0; j < currNodes.size(); j++ )
                        globalIDs[j] = currNodes[j].globalID();

                    // Get the libmesh element
                    const libMesh::Elem *currElemPtr =
                        libmeshElements.getElement( bnd1->globalID() );

                    // Get the DOF indicies for the matrix
                    for ( unsigned int i = 0; i < currNodes.size(); i++ )
                        d_dofManager->getDOFs( currNodes[i].globalID(), dofIndices[i] );

                    dofs.resize( currNodes.size() );
                    for ( size_t n = 0; n < dofIndices.size(); n++ )
                        dofs[n] = dofIndices[n][d_dofIds[nid][k]];

                    if ( d_isFluxGaussPtVector && myparams->d_variableFlux ) {
                        gpDOFManager->getDOFs( bnd1->globalID(), gpDofs );
                    }

                    d_fe->attach_quadrature_rule( d_qrule.get() );

                    d_phi = &( d_fe->get_phi() );
                    d_JxW = &( d_fe->get_JxW() );

                    d_fe->reinit( currElemPtr );

                    const auto &JxW          = *d_JxW;
                    const auto &phi          = *d_phi;
                    unsigned int numGaussPts = d_qrule->n_points();

                    std::vector<std::vector<double>> inputArgs(
                        elementInputVec.size(), std::vector<double>( currNodes.size() ) );
                    std::vector<std::vector<double>> inputArgsAtGpts(
                        elementInputVec.size(), std::vector<double>( numGaussPts ) );
                    std::vector<double> beta( numGaussPts, d_beta );
                    if ( d_robinPhysicsModel ) {
                        unsigned int startIdx = 0;
                        if ( d_isFluxGaussPtVector ) {
                            myparams->d_variableFlux->getValuesByGlobalID(
                                gpDofs.size(), &gpDofs[0], &inputArgsAtGpts[0][0] );
                            startIdx = 1;
                        }

                        for ( unsigned int m = startIdx; m < elementInputVec.size(); m++ ) {
                            elementInputVec[m]->getDOFManager()->getDOFs( globalIDs,
                                                                          dofsElementVec );
                            AMP_ASSERT( dofsElementVec.size() == dofIndices.size() );
                            elementInputVec[m]->getValuesByGlobalID(
                                dofsElementVec.size(), &dofsElementVec[0], &inputArgs[m][0] );
                            for ( size_t qp = 0; qp < currNodes.size(); qp++ ) {
                                for ( size_t n = 0; n < currNodes.size(); n++ ) {
                                    inputArgsAtGpts[m][qp] += phi[n][qp] * inputArgs[m][n];
                                }
                            }
                        }

                        std::vector<double> gamma( numGaussPts, d_gamma );
                        d_robinPhysicsModel->getConductance( beta, gamma, inputArgsAtGpts );
                    }

                    double temp;
                    for ( unsigned int qp = 0; qp < d_qrule->n_points(); qp++ ) {
                        for ( unsigned int j = 0; j < currNodes.size(); j++ ) {
                            for ( unsigned int i = 0; i < currNodes.size(); i++ ) {
                                temp = beta[qp] * ( JxW[qp] * phi[j][qp] * phi[i][qp] );
                                inputMatrix->addValueByGlobalID( dofs[j], dofs[i], temp );
                            } // end for i
                        }     // end for j
                    }         // end for qp

                } // end for bnd
            }     // end dof ids

        } // end for nid

        inputMatrix->makeConsistent();

    } // skip matrix
}
} // namespace AMP::Operator
