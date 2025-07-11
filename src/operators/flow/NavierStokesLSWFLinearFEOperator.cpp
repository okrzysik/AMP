#include "AMP/operators/flow/NavierStokesLSWFLinearFEOperator.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/node.h"

namespace AMP::Operator {

NavierStokesLSWFLinearFEOperator::NavierStokesLSWFLinearFEOperator(
    std::shared_ptr<const NavierStokesLinearFEOperatorParameters> params )
    : LinearFEOperator( params )
{
    AMP_INSIST( ( ( params.get() ) != nullptr ), "NULL parameter" );

    d_flowLSWFLinElem = std::dynamic_pointer_cast<NavierStokesLSWFLinearElement>( d_elemOp );

    AMP_INSIST( ( ( d_flowLSWFLinElem.get() ) != nullptr ),
                "d_elemOp is not of type NavierStokesLSWFLinearElement" );

    d_transportModel = params->d_transportModel;

    std::string varName = params->d_db->getString( "InputVariable" );
    d_inputVariable.reset( new AMP::LinearAlgebra::Variable( varName ) );
    d_outputVariable.reset( new AMP::LinearAlgebra::Variable( varName ) );

    /*
    std::vector<std::string> InternalVariableNames(4);
    InternalVariableNames[NavierStokes::PRESSURE]= "PRESSURE";
    InternalVariableNames[NavierStokes::VELOCITY]= "VELOCITY";
    InternalVariableNames[NavierStokes::PRINCIPALSTRESS]= "PRINCIPALSTRESS";
    InternalVariableNames[NavierStokes::SHEARSTRESS]= "SHEARSTRESS";

    AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
    auto activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");
    for(unsigned int i = 0; i < InternalVariableNames.size(); i++) {
        std::string varName = activeInpVar_db->getString(InternalVariableNames[i]);
        auto dummyVar = std::make_shared<AMP::LinearAlgebra::Variable>(varName);
        d_inpVariables->setVariable(i, dummyVar);
        d_outVariables->setVariable(i, dummyVar);
    }//end for i
    */

    bool isAttachedToNonlinearOperator =
        params->d_db->getWithDefault<bool>( "isAttachedToNonlinearOperator", false );

    if ( isAttachedToNonlinearOperator ) {
        bool isNonlinearOperatorInitialized =
            params->d_db->getWithDefault<bool>( "isNonlinearOperatorInitialized", false );
        if ( isNonlinearOperatorInitialized ) {
            reset( params );
        } else {
            AMP::LinearAlgebra::Vector::shared_ptr tmpInVec = AMP::LinearAlgebra::createVector(
                d_inDofMap, d_inputVariable, true, d_memory_location );
            AMP::LinearAlgebra::Vector::shared_ptr tmpOutVec = AMP::LinearAlgebra::createVector(
                d_outDofMap, d_outputVariable, true, d_memory_location );
            d_matrix = AMP::LinearAlgebra::createMatrix( tmpInVec, tmpOutVec );
        }
    } else {
        reset( params );
    }
}

void NavierStokesLSWFLinearFEOperator::preAssembly(
    std::shared_ptr<const OperatorParameters> oparams )
{

    auto params =
        std::dynamic_pointer_cast<const NavierStokesLinearFEOperatorParameters>( oparams );

    if ( params->d_frozenVec ) {
        d_inVec = mySubsetVector( params->d_frozenVec, d_inputVariable );
    }
    /*
          for(unsigned int i = 0; i < d_inpVariables->numVariables() ; i++) {
            std::shared_ptr<AMP::LinearAlgebra::Variable> var = d_inpVariables->getVariable(i);
            AMP::LinearAlgebra::Vector::shared_ptr vector = mySubsetVector(params->d_frozenVec[i],
       var);
            if((d_inVec[i].get() == NULL) and (vector.get() != NULL)) {
              d_inVec[i] = vector->clone();
            }
            if(d_inVec[i].get() != NULL) {
              if(params->d_frozenVec[i].get() != NULL) {
                d_inVec[i]->copyVector(vector);
                d_inVec[i]->makeConsistent(AMP::LinearAlgebra::ScatterType::CONSISTENT_SET);
              } else {
                d_inVec[i].reset();
              }
            }
          }
    */
    d_matrix->zero();
}

void NavierStokesLSWFLinearFEOperator::postAssembly()
{
    d_matrix->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
}

void NavierStokesLSWFLinearFEOperator::preElementOperation( const AMP::Mesh::MeshElement &elem )
{

    d_currNodes                     = elem.getElements( AMP::Mesh::GeomType::Vertex );
    unsigned int numNodesInCurrElem = d_currNodes.size();

    getDofIndicesForCurrentElement( NavierStokes::VELOCITY, d_type0DofIndices );
    //      getDofIndicesForCurrentElement(NavierStokes::PRESSURE, d_type1DofIndices);

    createCurrentLibMeshElement();

    std::vector<double> elementInputVectors;
    elementInputVectors.resize( 10 * numNodesInCurrElem );

    for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
        for ( unsigned int d = 0; d < 10; d++ ) {
            if ( d_inVec != nullptr ) {
                elementInputVectors[( 10 * r ) + d] =
                    d_inVec->getValueByGlobalID( d_type0DofIndices[r][d] );
            } else {
                elementInputVectors[( 10 * r ) + d] = 0.0;
            }
        } // end d
    }     // end r

    /*
          std::vector<std::vector<double> >
       elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

          elementInputVectors[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
          elementInputVectors[NavierStokes::PRESSURE].resize(numNodesInCurrElem);
          elementInputVectors[NavierStokes::PRINCIPALSTRESS].resize(3*numNodesInCurrElem);
          elementInputVectors[NavierStokes::SHEARSTRESS].resize(3*numNodesInCurrElem);

          for(unsigned int r = 0; r < numNodesInCurrElem ; r++) {
            for(unsigned int d = 0; d < 3; d++) {
              elementInputVectors[NavierStokes::VELOCITY][(3*r) + d]        =
       (d_inVec[NavierStokes::VELOCITY])->getValueByGlobalID( d_type0DofIndices[r][d] );
              elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*r) + d] =
       d_inVec[NavierStokes::PRINCIPALSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
              elementInputVectors[NavierStokes::SHEARSTRESS][(3*r) + d]     =
       d_inVec[NavierStokes::SHEARSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
            }

            elementInputVectors[NavierStokes::PRESSURE][r] =
       (d_inVec[NavierStokes::PRESSURE])->getValueByGlobalID(
       d_type1DofIndices[r][0] );

          }
    */

    unsigned int num_local_dofs = 10 * numNodesInCurrElem;
    d_elementStiffnessMatrix.resize( num_local_dofs );
    for ( unsigned int r = 0; r < num_local_dofs; r++ ) {
        d_elementStiffnessMatrix[r].resize( num_local_dofs );
        for ( unsigned int c = 0; c < num_local_dofs; c++ ) {
            d_elementStiffnessMatrix[r][c] = 0;
        }
    }

    d_flowLSWFLinElem->initializeForCurrentElement( d_currElemPtr, d_transportModel );
    d_flowLSWFLinElem->setElementVectors( elementInputVectors );
    d_flowLSWFLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
}

void NavierStokesLSWFLinearFEOperator::postElementOperation()
{
    PROFILE( "postElementOperation", 5 );
    for ( unsigned int r = 0; r < d_type0DofIndices.size(); r++ ) {
        for ( unsigned int dr = 0; dr < 10; dr++ ) {
            for ( unsigned int c = 0; c < d_type0DofIndices.size(); c++ ) {
                for ( unsigned int dc = 0; dc < 10; dc++ ) {
                    d_matrix->addValueByGlobalID(
                        d_type0DofIndices[r][dr],
                        d_type0DofIndices[c][dc],
                        d_elementStiffnessMatrix[( 10 * r ) + dr][( 10 * c ) + dc] );
                }
            }
        }
    }

    destroyCurrentLibMeshElement();
}

void NavierStokesLSWFLinearFEOperator::getDofIndicesForCurrentElement(
    int, std::vector<std::vector<size_t>> &dofIds )
{
    dofIds.resize( d_currNodes.size() );
    for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
        //        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds[j]);
        d_inDofMap->getDOFs( d_currNodes[j].globalID(), dofIds[j] );
    } // end of j
}

AMP::LinearAlgebra::Vector::shared_ptr NavierStokesLSWFLinearFEOperator::mySubsetVector(
    AMP::LinearAlgebra::Vector::shared_ptr vec, std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        auto meshSubsetVec = vec->select( meshSelector );
        return meshSubsetVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}
} // namespace AMP::Operator
