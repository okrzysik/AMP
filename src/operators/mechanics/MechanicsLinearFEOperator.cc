
#include "MechanicsLinearFEOperator.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"

namespace AMP {
namespace Operator {

MechanicsLinearFEOperator::MechanicsLinearFEOperator(
    const std::shared_ptr<MechanicsLinearFEOperatorParameters> &params )
    : LinearFEOperator( params )
{
    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameter" );
    d_useUpdatedLagrangian = ( params->d_db )->getWithDefault( "USE_UPDATED_LAGRANGIAN", false );
    if ( d_useUpdatedLagrangian ) {
        d_mechLinULElem =
            std::dynamic_pointer_cast<MechanicsLinearUpdatedLagrangianElement>( d_elemOp );
        AMP_INSIST( ( ( d_mechLinULElem.get() ) != nullptr ),
                    "d_elemOp is not of type MechanicsLinearUpdatedLagrangianElement" );
    } else {
        d_mechLinElem = std::dynamic_pointer_cast<MechanicsLinearElement>( d_elemOp );
        AMP_INSIST( ( ( d_mechLinElem.get() ) != nullptr ),
                    "d_elemOp is not of type MechanicsLinearElement" );
    }
    d_materialModel = params->d_materialModel;

    AMP_INSIST( params->d_db->keyExists( "InputVariable" ), "key not found" );
    std::string inpVarName = params->d_db->getString( "InputVariable" );
    d_inputVariable.reset( new AMP::LinearAlgebra::Variable( inpVarName ) );

    AMP_INSIST( params->d_db->keyExists( "OutputVariable" ), "key not found" );
    std::string outVarName = params->d_db->getString( "OutputVariable" );
    d_outputVariable.reset( new AMP::LinearAlgebra::Variable( outVarName ) );

    if ( d_useUpdatedLagrangian ) {
        d_refXYZ = AMP::LinearAlgebra::createVector( d_inDofMap, d_inputVariable, true );
        d_refXYZ->zero();

        AMP::Mesh::MeshIterator el     = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
        AMP::Mesh::MeshIterator end_el = el.end();

        for ( ; el != end_el; ++el ) {
            d_currNodes                     = el->getElements( AMP::Mesh::GeomType::Vertex );
            unsigned int numNodesInCurrElem = d_currNodes.size();

            getDofIndicesForCurrentElement();

            createCurrentLibMeshElement();

            std::vector<double> elementRefXYZ;
            elementRefXYZ.resize( 3 * numNodesInCurrElem );

            d_mechLinULElem->initializeForCurrentElement( d_currElemPtr, d_materialModel );
            d_mechLinULElem->initializeReferenceXYZ( elementRefXYZ );

            for ( unsigned int j = 0; j < numNodesInCurrElem; j++ ) {
                for ( unsigned int i = 0; i < 3; i++ ) {
                    d_refXYZ->setValueByGlobalID( d_dofIndices[j][i],
                                                  elementRefXYZ[( 3 * j ) + i] );
                } // end of i
            }     // end of j

            destroyCurrentLibMeshElement();
        } // end of el

        d_refXYZ->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    } // end of UpdatedLagrangian condition.

    bool isAttachedToNonlinearOperator =
        params->d_db->getWithDefault( "isAttachedToNonlinearOperator", false );

    if ( isAttachedToNonlinearOperator ) {
        bool isNonlinearOperatorInitialized =
            params->d_db->getWithDefault( "isNonlinearOperatorInitialized", false );
        if ( isNonlinearOperatorInitialized ) {
            reset( params );
        } else {
            AMP::LinearAlgebra::Vector::shared_ptr tmpInVec =
                AMP::LinearAlgebra::createVector( d_inDofMap, d_inputVariable, true );
            AMP::LinearAlgebra::Vector::shared_ptr tmpOutVec =
                AMP::LinearAlgebra::createVector( d_outDofMap, d_outputVariable, true );
            d_matrix = AMP::LinearAlgebra::createMatrix( tmpInVec, tmpOutVec );
        }
    } else {
        reset( params );
    }
}

void MechanicsLinearFEOperator::preAssembly( const std::shared_ptr<OperatorParameters> &oparams )
{
    if ( d_useUpdatedLagrangian ) {
        std::shared_ptr<MechanicsLinearFEOperatorParameters> params =
            std::dynamic_pointer_cast<MechanicsLinearFEOperatorParameters>( oparams );
        AMP_INSIST( ( params != nullptr ), "NULL params" );

        if ( ( d_dispVec == nullptr ) and ( params->d_dispVec != nullptr ) ) {
            d_dispVec = ( params->d_dispVec )->cloneVector();
        }
        if ( d_dispVec != nullptr ) {
            if ( params->d_dispVec != nullptr ) {
                d_dispVec->copyVector( params->d_dispVec );
                d_dispVec->makeConsistent(
                    AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
            } else {
                d_dispVec.reset();
            }
        }
    }
    d_matrix->zero();
    d_materialModel->preLinearAssembly();
}

void MechanicsLinearFEOperator::preElementOperation( const AMP::Mesh::MeshElement &elem )
{
    d_currNodes                     = elem.getElements( AMP::Mesh::GeomType::Vertex );
    unsigned int numNodesInCurrElem = d_currNodes.size();

    getDofIndicesForCurrentElement();

    createCurrentLibMeshElement();

    unsigned int num_local_dofs = ( 3 * numNodesInCurrElem );
    d_elementStiffnessMatrix.resize( num_local_dofs );
    for ( unsigned int r = 0; r < num_local_dofs; r++ ) {
        d_elementStiffnessMatrix[r].resize( num_local_dofs );
        for ( unsigned int c = 0; c < num_local_dofs; c++ ) {
            d_elementStiffnessMatrix[r][c] = 0.0;
        } // end for c
    }     // end for r

    if ( d_useUpdatedLagrangian ) {
        std::vector<std::vector<double>> elementInputVectors(
            Mechanics::TOTAL_NUMBER_OF_VARIABLES );
        std::vector<double> elementRefXYZ;
        elementInputVectors[Mechanics::DISPLACEMENT].resize( num_local_dofs );
        elementRefXYZ.resize( num_local_dofs );
        for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
            for ( unsigned int d = 0; d < 3; d++ ) {
                if ( d_dispVec != nullptr ) {
                    elementInputVectors[Mechanics::DISPLACEMENT][( 3 * r ) + d] =
                        d_dispVec->getValueByGlobalID( d_dofIndices[r][d] );
                } else {
                    elementInputVectors[Mechanics::DISPLACEMENT][( 3 * r ) + d] = 0.0;
                }
                elementRefXYZ[( 3 * r ) + d] = d_refXYZ->getValueByGlobalID( d_dofIndices[r][d] );
            } // end for d
        }     // end for r
        d_mechLinULElem->initializeForCurrentElement( d_currElemPtr, d_materialModel );
        d_mechLinULElem->setElementVectors( elementInputVectors );
        d_mechLinULElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
        d_mechLinULElem->assignReferenceXYZ( elementRefXYZ );
    } else {
        d_mechLinElem->initializeForCurrentElement( d_currElemPtr, d_materialModel );
        d_mechLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
    }
}

void MechanicsLinearFEOperator::postElementOperation()
{
    unsigned int numNodesInCurrElem = d_dofIndices.size();
    for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
        for ( unsigned int dr = 0; dr < 3; dr++ ) {
            for ( unsigned int c = 0; c < numNodesInCurrElem; c++ ) {
                for ( unsigned int dc = 0; dc < 3; dc++ ) {
                    d_matrix->addValueByGlobalID(
                        d_dofIndices[r][dr],
                        d_dofIndices[c][dc],
                        d_elementStiffnessMatrix[( 3 * r ) + dr][( 3 * c ) + dc] );
                } // end for dc
            }     // end for c
        }         // end for dr
    }             // end for r
    destroyCurrentLibMeshElement();
}

void MechanicsLinearFEOperator::postAssembly()
{
    d_materialModel->postLinearAssembly();
    d_matrix->makeConsistent();
}

void MechanicsLinearFEOperator::printStressAndStrain( AMP::LinearAlgebra::Vector::shared_ptr disp,
                                                      const std::string &fname )
{
    FILE *fp = fopen( fname.c_str(), "w" );
    fprintf( fp,
             "x, y, z, Stresses(11, 22, 33, 23, 13, 12), Strains(11, 22, 33, 23, 13, 12) \n\n" );

    disp->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    d_materialModel->preLinearAssembly();

    AMP::Mesh::MeshIterator el     = d_Mesh->getIterator( AMP::Mesh::GeomType::Volume, 0 );
    AMP::Mesh::MeshIterator end_el = el.end();

    for ( ; el != end_el; ++el ) {
        d_currNodes                     = el->getElements( AMP::Mesh::GeomType::Vertex );
        unsigned int numNodesInCurrElem = d_currNodes.size();

        getDofIndicesForCurrentElement();

        createCurrentLibMeshElement();

        std::vector<double> elementInputVector( 3 * numNodesInCurrElem );
        for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
            for ( unsigned int d = 0; d < 3; d++ ) {
                elementInputVector[( 3 * r ) + d] = disp->getValueByGlobalID( d_dofIndices[r][d] );
            } // end for d
        }     // end for r

        d_mechLinElem->initializeForCurrentElement( d_currElemPtr, d_materialModel );
        d_mechLinElem->printStressAndStrain( fp, elementInputVector );

        destroyCurrentLibMeshElement();
    } // end for el

    d_materialModel->postLinearAssembly();
    fprintf( fp, "\n\n" );
    fclose( fp );
}

void MechanicsLinearFEOperator::getDofIndicesForCurrentElement()
{
    d_dofIndices.resize( d_currNodes.size() );
    for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
        d_inDofMap->getDOFs( d_currNodes[j].globalID(), d_dofIndices[j] );
    } // end of j
}
} // namespace Operator
} // namespace AMP
