
#include "NavierStokesLSWFFEOperator.h"
#include "NavierStokesLSWFFEOperatorParameters.h"
#include "NavierStokesLinearFEOperatorParameters.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "vectors/VectorBuilder.h"
#include "cell_hex8.h"
#include "node.h"

namespace AMP {
  namespace Operator {

    NavierStokesLSWFFEOperator :: NavierStokesLSWFFEOperator (
        const boost::shared_ptr<NavierStokesLSWFFEOperatorParameters> & params)
      : NonlinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
        AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

        d_nsLSWFElem = boost::dynamic_pointer_cast<NavierStokesLSWFElement>(d_elemOp);

        AMP_INSIST( ((d_nsLSWFElem.get()) != NULL), "d_elemOp is not of type NavierStokesLSWFElement" );

        d_transportModel = params->d_transportModel;

//        d_isActive.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES, true);
//        d_isFrozen.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES, false);
//        d_inVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

//        d_isActive[NavierStokes::TEMPERATURE] = false;
//        d_isFrozen[NavierStokes::TEMPERATURE] = false;

//        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          d_dofMap = (params->d_dofMap);
//        }//end for i

/*
        d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
          d_inpVariables->add(dummyVar);
        }//end for i
        AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
        boost::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

        std::vector<std::string> InternalVariableNames(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
        InternalVariableNames[NavierStokes::PRESSURE]= "PRESSURE";
        InternalVariableNames[NavierStokes::VELOCITY]= "VELOCITY";
        InternalVariableNames[NavierStokes::PRINCIPALSTRESS]= "PRINCIPALSTRESS";
        InternalVariableNames[NavierStokes::SHEARSTRESS]= "SHEARSTRESS";
        InternalVariableNames[NavierStokes::TEMPERATURE]= "TEMPERATURE";

        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          if(d_isActive[i]) {
*/          
            std::string varName = params->d_db->getString("InputVariable");
            d_inpVariables.reset(new AMP::LinearAlgebra::Variable(varName) );
            d_outVariables.reset(new AMP::LinearAlgebra::Variable(varName) );
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

    void NavierStokesLSWFFEOperator :: preAssembly(const AMP::LinearAlgebra::Vector::shared_ptr &u, 
        boost::shared_ptr< AMP::LinearAlgebra::Vector >  &r) {
      AMP_INSIST( (u != NULL), "NULL Input Vector" );
/*
      for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
        if(d_isActive[i]) {
          if(!(d_isFrozen[i])) {
            AMP::LinearAlgebra::Variable::shared_ptr var = d_inpVariables->getVariable(i); 
            AMP::LinearAlgebra::Vector::shared_ptr vector = mySubsetVector(u, d_inpVariables);
*/          
            d_inVec = mySubsetVector(u, d_inpVariables);
/*            
            setVector(i, vector);
          }
        }
      }//end for i
*/
      d_outVec = mySubsetVector(r, d_outVariables);
      d_outVec->zero();

    }

    void NavierStokesLSWFFEOperator :: postAssembly() {
      d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    void NavierStokesLSWFFEOperator :: preElementOperation( const AMP::Mesh::MeshElement & elem ) {
      d_currNodes = elem.getElements(AMP::Mesh::Vertex);
      unsigned int numNodesInCurrElem = d_currNodes.size();

      getDofIndicesForCurrentElement(NavierStokes::VELOCITY, d_type0DofIndices);
//      getDofIndicesForCurrentElement(NavierStokes::PRESSURE, d_type1DofIndices);

      std::vector<double> elementInputVectors;
      elementInputVectors.resize(10*numNodesInCurrElem);
/*      
      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
      elementInputVectors[NavierStokes::PRESSURE].resize(numNodesInCurrElem);
      elementInputVectors[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
      elementInputVectors[NavierStokes::PRINCIPALSTRESS].resize(3*numNodesInCurrElem);
      elementInputVectors[NavierStokes::SHEARSTRESS].resize(3*numNodesInCurrElem);
      if(d_isActive[NavierStokes::TEMPERATURE]) {
        elementInputVectors[NavierStokes::TEMPERATURE].resize(numNodesInCurrElem);
      }

      for(unsigned int r = 0; r < numNodesInCurrElem; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d]        = d_inVec[NavierStokes::VELOCITY]->getValueByGlobalID( d_type0DofIndices[r][d] );
          elementInputVectors[NavierStokes::PRINCIPALSTRESS][(3*r) + d] = d_inVec[NavierStokes::PRINCIPALSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
          elementInputVectors[NavierStokes::SHEARSTRESS][(3*r) + d]     = d_inVec[NavierStokes::SHEARSTRESS]->getValueByGlobalID( d_type0DofIndices[r][d] );
        }//end d
        elementInputVectors[NavierStokes::PRESSURE][r]                  = d_inVec[NavierStokes::PRESSURE]->getValueByGlobalID( d_type1DofIndices[r][0] );
        if(d_isActive[NavierStokes::TEMPERATURE]) {
          elementInputVectors[NavierStokes::TEMPERATURE][r]             = d_inVec[NavierStokes::TEMPERATURE]->getValueByGlobalID( d_type1DofIndices[r][0] );
        }
      }//end r
*/
 
      for(unsigned int r = 0; r < numNodesInCurrElem; r++) {
        for(unsigned int d = 0; d < 10; d++) {
          elementInputVectors[(10*r) + d]       = d_inVec->getValueByGlobalID( d_type0DofIndices[r][d]   );
        }//end d
      }//end r     

      d_elementOutputVector.resize(10*numNodesInCurrElem);
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
      d_nsLSWFElem->setElementVectors( elementInputVectors , d_elementOutputVector );
    }

    void NavierStokesLSWFFEOperator :: postElementOperation() {

/*      
      AMP::LinearAlgebra::Variable::shared_ptr velocityVar        = d_outVariables->getVariable(NavierStokes::VELOCITY); 
      AMP::LinearAlgebra::Variable::shared_ptr pressureVar        = d_outVariables->getVariable(NavierStokes::PRESSURE); 
      AMP::LinearAlgebra::Variable::shared_ptr principalStressVar = d_outVariables->getVariable(NavierStokes::PRINCIPALSTRESS); 
      AMP::LinearAlgebra::Variable::shared_ptr shearStressVar     = d_outVariables->getVariable(NavierStokes::SHEARSTRESS); 
      AMP::LinearAlgebra::Vector::shared_ptr d_velOutVec = mySubsetVector(d_outVec, velocityVar );
      AMP::LinearAlgebra::Vector::shared_ptr d_preOutVec = mySubsetVector(d_outVec, pressureVar );
      AMP::LinearAlgebra::Vector::shared_ptr d_pstOutVec = mySubsetVector(d_outVec, principalStressVar );
      AMP::LinearAlgebra::Vector::shared_ptr d_sstOutVec = mySubsetVector(d_outVec, shearStressVar );
      
      for(unsigned int r = 0; r < d_type0DofIndices.size(); r++) {
        for(unsigned int d = 0; d < 3; d++) {
          d_velOutVec->addValueByGlobalID( d_type0DofIndices[r][d], d_elementOutputVector[NavierStokes::VELOCITY][(3*r) + d] );
          d_pstOutVec->addValueByGlobalID( d_type0DofIndices[r][d], d_elementOutputVector[NavierStokes::PRINCIPALSTRESS][(3*r) + d] );
          d_sstOutVec->addValueByGlobalID( d_type0DofIndices[r][d], d_elementOutputVector[NavierStokes::SHEARSTRESS][(3*r) + d] );
        }//end for d
      }//end for r
      for(unsigned int r = 0; r < d_type1DofIndices.size(); r++) {
          d_preOutVec->addValueByGlobalID( d_type1DofIndices[r][0], d_elementOutputVector[NavierStokes::PRESSURE][r] );
      }
*/

      for(unsigned int r = 0; r < d_type0DofIndices.size(); r++) {
        AMP_ASSERT(d_type0DofIndices[r].size() == 10);
        for(unsigned int d = 0; d < 10; d++) {
          d_outVec->addValueByGlobalID( d_type0DofIndices[r][d], d_elementOutputVector[(10*r) + d] );
        }//end for d
      }//end for r

    }

    void NavierStokesLSWFFEOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
       // DO Nothing
    }

    boost::shared_ptr<OperatorParameters> NavierStokesLSWFFEOperator ::
      getJacobianParameters(const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) {

        // set up a database for the linear operator params
        boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
        tmp_db->putBool("reset_reuses_matrix", true);
        tmp_db->putBool("isAttachedToNonlinearOperator", true);

        // create the linear operator params
        boost::shared_ptr<NavierStokesLinearFEOperatorParameters> outParams(new
            NavierStokesLinearFEOperatorParameters(tmp_db));
/*
        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          if(d_isActive[i]) {
              AMP::LinearAlgebra::Variable::shared_ptr var = d_inpVariables->getVariable(i); 
              AMP::LinearAlgebra::Vector::shared_ptr vector = mySubsetVector(u, var);
              outParams->d_frozenVec[i] = vector ;
              (outParams->d_frozenVec[i])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
          }
        }//end for i
*/
        outParams->d_frozenVec = mySubsetVector(u, d_inpVariables);
        outParams->d_frozenVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        d_outVec.reset();

        return outParams;
      }

    void NavierStokesLSWFFEOperator :: getDofIndicesForCurrentElement(int varId, std::vector<std::vector<size_t> > & dofIds) {
      dofIds.resize(d_currNodes.size());
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
//        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds[j]);
        d_dofMap->getDOFs(d_currNodes[j].globalID(), dofIds[j]);
      } // end of j
    }

  }
}//end namespace

