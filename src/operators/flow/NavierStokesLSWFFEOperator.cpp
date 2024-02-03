#include "NavierStokesLSWFFEOperator.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#include "NavierStokesLSWFFEOperatorParameters.h"
#include "NavierStokesLinearFEOperatorParameters.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/node.h"

namespace AMP::Operator {

NavierStokesLSWFFEOperator::NavierStokesLSWFFEOperator(
    std::shared_ptr<const NavierStokesLSWFFEOperatorParameters> params )
    : NonlinearFEOperator( params )
{
    AMP_INSIST( params, "NULL parameter!" );
    AMP_INSIST( params->d_db, "NULL database!" );

    d_nsLSWFElem = std::dynamic_pointer_cast<NavierStokesLSWFElement>( d_elemOp );

    AMP_INSIST( ( ( d_nsLSWFElem.get() ) != nullptr ),
                "d_elemOp is not of type NavierStokesLSWFElement" );

    d_transportModel = params->d_transportModel;

    // d_isActive.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES, true);
    // d_isFrozen.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES, false);
    // d_inVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
    // d_isActive[NavierStokes::TEMPERATURE] = false;
    // d_isFrozen[NavierStokes::TEMPERATURE] = false;

    // for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
    d_dofMap = ( params->d_dofMap );
    // }//end for i

    /*
    d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
    for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
      std::shared_ptr<AMP::LinearAlgebra::Variable> dummyVar;
      d_inpVariables->add(dummyVar);
    }//end for i
    AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
    auto activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");
    std::vector<std::string> InternalVariableNames(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
    InternalVariableNames[NavierStokes::PRESSURE]= "PRESSURE";
    InternalVariableNames[NavierStokes::VELOCITY]= "VELOCITY";
    InternalVariableNames[NavierStokes::PRINCIPALSTRESS]= "PRINCIPALSTRESS";
    InternalVariableNames[NavierStokes::SHEARSTRESS]= "SHEARSTRESS";
    InternalVariableNames[NavierStokes::TEMPERATURE]= "TEMPERATURE";

    for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
      if(d_isActive[i]) {
    */
    std::string varName = params->d_db->getString( "InputVariable" );
    d_inpVariables.reset( new AMP::LinearAlgebra::Variable( varName ) );
    d_outVariables.reset( new AMP::LinearAlgebra::Variable( varName ) );
    /*
        d_inpVariables->setVariable(i, dummyVar);
        d_outVariables->setVariable(i, dummyVar);
        if(d_isFrozen[i]) {
          if( params->d_frozenVec[i] != NULL ) {
            setVector(i, params->d_frozenVec[i]);
          }
        }
      }
    }//end for i
    */
}

void NavierStokesLSWFFEOperator::preAssembly( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                              AMP::LinearAlgebra::Vector::shared_ptr r )
{
    AMP_INSIST( ( u != nullptr ), "NULL Input Vector" );
    /*
          for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
            if(d_isActive[i]) {
              if(!(d_isFrozen[i])) {
                std::shared_ptr<AMP::LinearAlgebra::Variable> var = d_inpVariables->getVariable(i);
                AMP::LinearAlgebra::Vector::shared_ptr vector = mySubsetVector(u, d_inpVariables);
    */
    d_inVec = mySubsetVector( u, d_inpVariables );
    /*
                setVector(i, vector);
              }
            }
          }//end for i
    */
    d_outVec = mySubsetVector( r, d_outVariables );
    d_outVec->zero();
}

void NavierStokesLSWFFEOperator::postAssembly()
{
    d_outVec->makeConsistent( AMP::LinearAlgebra::ScatterType::CONSISTENT_ADD );
}

void NavierStokesLSWFFEOperator::preElementOperation( const AMP::Mesh::MeshElement &elem )
{
    d_currNodes                     = elem.getElements( AMP::Mesh::GeomType::Vertex );
    unsigned int numNodesInCurrElem = d_currNodes.size();

    getDofIndicesForCurrentElement( NavierStokes::VELOCITY, d_type0DofIndices );
    //      getDofIndicesForCurrentElement(NavierStokes::PRESSURE, d_type1DofIndices);

    std::vector<double> elementInputVectors;
    elementInputVectors.resize( 10 * numNodesInCurrElem );
    /*
    auto elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
    elementInputVectors[NavierStokes::PRESSURE].resize(numNodesInCurrElem);
    elementInputVectors[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
    elementInputVectors[NavierStokes::PRINCIPALSTRESS].resize(3*numNodesInCurrElem);
    elementInputVectors[NavierStokes::SHEARSTRESS].resize(3*numNodesInCurrElem);
    if(d_isActive[NavierStokes::TEMPERATURE]) {
        elementInputVectors[NavierStokes::TEMPERATURE].resize(numNodesInCurrElem);
    }

    for(unsigned int r = 0; r < numNodesInCurrElem; r++) {
       for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d]        =
             d_inVec[NavierStokes::VELOCITY]->getValueByGlobalID( d_type0DofIndices[r][d] );
          elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*r) + d] =
             d_inVec[NavierStokes::PRINCIPALSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
          elementInputVectors[NavierStokes::SHEARSTRESS][(3*r) + d]     =
             d_inVec[NavierStokes::SHEARSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
       }//end d
       elementInputVectors[NavierStokes::PRESSURE][r]                  =
           d_inVec[NavierStokes::PRESSURE]->getValueByGlobalID( d_type1DofIndices[r][0] );
       if(d_isActive[NavierStokes::TEMPERATURE]) {
           elementInputVectors[NavierStokes::TEMPERATURE][r]             =
              d_inVec[NavierStokes::TEMPERATURE]->getValueByGlobalID( d_type1DofIndices[r][0] );
       }
    }//end r
    */

    for ( unsigned int r = 0; r < numNodesInCurrElem; r++ ) {
        for ( unsigned int d = 0; d < 10; d++ ) {
            elementInputVectors[( 10 * r ) + d] =
                d_inVec->getValueByGlobalID( d_type0DofIndices[r][d] );
        } // end d
    }     // end r

    d_elementOutputVector.resize( 10 * numNodesInCurrElem );
    /*
         d_elementOutputVector.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
         d_elementOutputVector[NavierStokes::PRESSURE].resize(numNodesInCurrElem);
         d_elementOutputVector[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
         d_elementOutputVector[NavierStokes::PRINCIPALSTRESS].resize(3*numNodesInCurrElem);
         d_elementOutputVector[NavierStokes::SHEARSTRESS].resize(3*numNodesInCurrElem);
         if(d_isActive[4]) {
           d_elementOutputVector[4].resize(numNodesInCurrElem);
         }
   */

    d_nsLSWFElem->initializeForCurrentElement( d_currElemPtrs[d_currElemIdx], d_transportModel );
    d_nsLSWFElem->setElementVectors( elementInputVectors, d_elementOutputVector );
}

void NavierStokesLSWFFEOperator::postElementOperation()
{

    /*
    auto velocityVar = d_outVariables->getVariable(NavierStokes::VELOCITY);
    auto pressureVar = d_outVariables->getVariable(NavierStokes::PRESSURE);
    auto principalStressVar = d_outVariables->getVariable(NavierStokes::PRINCIPALSTRESS);
    auto shearStressVar = d_outVariables->getVariable(NavierStokes::SHEARSTRESS);
    auto d_velOutVec = mySubsetVector(d_outVec, velocityVar );
    auto d_preOutVec = mySubsetVector(d_outVec, pressureVar );
    auto d_pstOutVec = mySubsetVector(d_outVec, principalStressVar );
    auto d_sstOutVec = mySubsetVector(d_outVec, shearStressVar );
    for(unsigned int r = 0; r < d_type0DofIndices.size(); r++) {
        for(unsigned int d = 0; d < 3; d++) {
           d_velOutVec->addValueByGlobalID( d_type0DofIndices[r][d],
              d_elementOutputVector[NavierStokes::VELOCITY][(3*r) + d] );
           d_pstOutVec->addValueByGlobalID( d_type0DofIndices[r][d],
              d_elementOutputVector[NavierStokes::PRINCIPALSTRESS][(3*r) + d] );
           d_sstOutVec->addValueByGlobalID( d_type0DofIndices[r][d],
              d_elementOutputVector[NavierStokes::SHEARSTRESS][(3*r) + d] );
        }//end for d
    }//end for r
    for(unsigned int r = 0; r < d_type1DofIndices.size(); r++) {
        d_preOutVec->addValueByGlobalID( d_type1DofIndices[r][0],
           d_elementOutputVector[NavierStokes::PRESSURE][r] );
    }
    */

    for ( unsigned int r = 0; r < d_type0DofIndices.size(); r++ ) {
        AMP_ASSERT( d_type0DofIndices[r].size() == 10 );
        for ( unsigned int d = 0; d < 10; d++ ) {
            d_outVec->addValuesByGlobalID(
                1, &d_type0DofIndices[r][d], &d_elementOutputVector[( 10 * r ) + d] );
        } // end for d
    }     // end for r
}

void NavierStokesLSWFFEOperator::reset( std::shared_ptr<const OperatorParameters> )
{
    // DO Nothing
}

std::shared_ptr<OperatorParameters> NavierStokesLSWFFEOperator::getJacobianParameters(
    AMP::LinearAlgebra::Vector::const_shared_ptr u_in )
{

    auto u = std::const_pointer_cast<AMP::LinearAlgebra::Vector>( u_in );

    // set up a database for the linear operator params
    auto tmp_db = std::make_shared<AMP::Database>( "Dummy" );
    tmp_db->putScalar( "reset_reuses_matrix", true );
    tmp_db->putScalar( "isAttachedToNonlinearOperator", true );

    // create the linear operator params
    auto outParams         = std::make_shared<NavierStokesLinearFEOperatorParameters>( tmp_db );
    outParams->d_frozenVec = mySubsetVector( u, d_inpVariables );
    outParams->d_frozenVec->makeConsistent(
        AMP::LinearAlgebra::ScatterType::CONSISTENT_SET );
    d_outVec.reset();

    return outParams;
}

void NavierStokesLSWFFEOperator::getDofIndicesForCurrentElement(
    int, std::vector<std::vector<size_t>> &dofIds )
{
    dofIds.resize( d_currNodes.size() );
    for ( unsigned int j = 0; j < d_currNodes.size(); j++ ) {
        d_dofMap->getDOFs( d_currNodes[j].globalID(), dofIds[j] );
    }
}

AMP::LinearAlgebra::Vector::shared_ptr
NavierStokesLSWFFEOperator::mySubsetVector( AMP::LinearAlgebra::Vector::shared_ptr vec,
                                            std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        AMP::LinearAlgebra::Vector::shared_ptr meshSubsetVec =
            vec->select( meshSelector, var->getName() );
        return meshSubsetVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}

AMP::LinearAlgebra::Vector::const_shared_ptr
NavierStokesLSWFFEOperator::mySubsetVector( AMP::LinearAlgebra::Vector::const_shared_ptr vec,
                                            std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
    if ( d_Mesh ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector( d_Mesh );
        AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetVec =
            vec->select( meshSelector, var->getName() );
        return meshSubsetVec->subsetVectorForVariable( var );
    } else {
        return vec->subsetVectorForVariable( var );
    }
}
} // namespace AMP::Operator
