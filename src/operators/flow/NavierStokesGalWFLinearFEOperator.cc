
#include "operators/flow/NavierStokesGalWFLinearFEOperator.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/VectorBuilder.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/node.h"

namespace AMP {
  namespace Operator {
/*
    NavierStokesGalWFLinearFEOperator :: NavierStokesGalWFLinearFEOperator (
        const AMP::shared_ptr<NavierStokesLinearFEOperatorParameters> & params)
      : LinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

        d_flowGalWFLinElem = AMP::dynamic_pointer_cast<NavierStokesGalWFLinearElement>(d_elemOp);

        AMP_INSIST( ((d_flowGalWFLinElem.get()) != NULL), "d_elemOp is not of type NavierStokesGalWFLinearElement" );


        d_transportModel = params->d_transportModel ;

        d_inpVariable.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
        d_outVariable.reset(new AMP::LinearAlgebra::MultiVariable("myOutVar"));
        for(unsigned int i = 0; i < 2; i++) {
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
          d_inpVariable->add(dummyVar);
          d_outVariable->add(dummyVar);
        }//end for i

        AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
        AMP::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

        AMP_INSIST(activeInpVar_db->keyExists("VELOCITY"), "VELOCITY must be active");
        AMP_INSIST(activeInpVar_db->keyExists("PRESSURE"), "VELOCITY must be active");

        std::string tempVarName = activeInpVar_db->getString("VELOCITY");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::Variable(tempVarName) ); 
        d_inpVariable->setVariable(NavierStokes::VELOCITY, tempVar);
        d_outVariable->setVariable(NavierStokes::VELOCITY, tempVar);

        tempVarName = activeInpVar_db->getString("PRESSURE");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar2(new AMP::LinearAlgebra::Variable(tempVarName) ); 
        d_inpVariable->setVariable(NavierStokes::PRESSURE, tempVar2);
        d_outVariable->setVariable(NavierStokes::PRESSURE, tempVar2);

        bool isAttachedToNonlinearOperator = params->d_db->getBoolWithDefault("isAttachedToNonlinearOperator", false);

        if(isAttachedToNonlinearOperator) {
          bool isNonlinearOperatorInitialized = params->d_db->getBoolWithDefault("isNonlinearOperatorInitialized", false);
          if(isNonlinearOperatorInitialized) {
            reset(params);
          } else {
            AMP::LinearAlgebra::Vector::shared_ptr tmpInVec = AMP::LinearAlgebra::createVector(d_inDofMap, d_inpVariable, true);
            AMP::LinearAlgebra::Vector::shared_ptr tmpOutVec = AMP::LinearAlgebra::createVector(d_outDofMap, d_outVariable, true);
            d_matrix = AMP::LinearAlgebra::createMatrix(tmpInVec, tmpOutVec);
          }
        } else {
          reset(params);
        }
      }

    void NavierStokesGalWFLinearFEOperator :: preAssembly(const AMP::shared_ptr<OperatorParameters>& oparams) 
    {

      AMP::shared_ptr<NavierStokesLinearFEOperatorParameters> params = AMP::dynamic_pointer_cast<NavierStokesLinearFEOperatorParameters>(oparams);

      if((d_inVec[NavierStokes::VELOCITY].get() == NULL) and (params->d_frozenVec[NavierStokes::VELOCITY].get() != NULL)) {
        d_inVec[NavierStokes::VELOCITY] = params->d_frozenVec[NavierStokes::VELOCITY]->cloneVector();
      }
      if(d_inVec[NavierStokes::VELOCITY].get() != NULL) {
        if(params->d_frozenVec[NavierStokes::VELOCITY].get() != NULL) {
          d_inVec[NavierStokes::VELOCITY]->copyVector(params->d_frozenVec[NavierStokes::VELOCITY]);
          d_inVec[NavierStokes::VELOCITY]->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
        } else {
          d_inVec[NavierStokes::VELOCITY].reset();
        }
      }

      if((d_inVec[NavierStokes::PRESSURE].get() == NULL) and (params->d_frozenVec[NavierStokes::PRESSURE].get() != NULL)) {
        d_inVec[NavierStokes::PRESSURE] = params->d_frozenVec[NavierStokes::PRESSURE]->cloneVector();
      }
      if(d_inVec[NavierStokes::PRESSURE].get() != NULL) {
        if(params->d_frozenVec[NavierStokes::PRESSURE].get() != NULL) {
          d_inVec[NavierStokes::PRESSURE]->copyVector(params->d_frozenVec[NavierStokes::PRESSURE]);
          d_inVec[NavierStokes::PRESSURE]->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
        } else {
          d_inVec[NavierStokes::PRESSURE].reset();
        }
      }

      if(d_inVec[NavierStokes::TEMPERATURE].get() == NULL and params->d_frozenVec[NavierStokes::TEMPERATURE].get() != NULL) {
        d_inVec[NavierStokes::TEMPERATURE] = params->d_frozenVec[NavierStokes::TEMPERATURE]->cloneVector();
      }
      if(d_inVec[NavierStokes::TEMPERATURE].get() != NULL) {
        if(params->d_frozenVec[NavierStokes::TEMPERATURE].get() != NULL) {
          d_inVec[NavierStokes::TEMPERATURE]->copyVector(params->d_frozenVec[NavierStokes::TEMPERATURE]);
          d_inVec[NavierStokes::TEMPERATURE]->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);
        } else {
          d_inVec[NavierStokes::TEMPERATURE].reset();
        }
      }


      d_matrix->zero();

    }

    void NavierStokesGalWFLinearFEOperator :: postAssembly()
    {

      d_matrix->makeConsistent ();
    }

    void NavierStokesGalWFLinearFEOperator :: preElementOperation( const AMP::Mesh::MeshElement & elem ) {

      d_currNodes = elem.getElements(AMP::Mesh::Vertex);
      unsigned int numNodesInCurrElem = d_currNodes.size();

      gettype0DofIndicesForCurrentElement(NavierStokes::VELOCITY, d_type0DofIndices);
      gettype1DofIndicesForCurrentElement(NavierStokes::PRESSURE, d_type1DofIndices);

      createHex27LibMeshElement();

      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[NavierStokes::VELOCITY].resize(3*numNodesInCurrElem);
      elementInputVectors[NavierStokes::PRESSURE].resize(numNodesInCurrElem);

      if(d_inVec[NavierStokes::TEMPERATURE].get() != NULL) {
        elementInputVectors[NavierStokes::TEMPERATURE].resize(numNodesInCurrElem);
      }

      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->getValueByGlobalID( d_type0DofIndices[r][d] );
        }

        elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->getValueByGlobalID( d_type1DofIndices[r] );

        if(d_inVec[NavierStokes::TEMPERATURE].get() != NULL) {
          elementInputVectors[NavierStokes::TEMPERATURE][r] = (d_inVec[NavierStokes::TEMPERATURE])->getValueByGlobalID( d_type1DofIndices[r] );
        }
      }


      unsigned int num_local_dofs = 4*numNodesInCurrElem;
      d_elementStiffnessMatrix.resize(num_local_dofs);
      for(unsigned int r = 0; r < num_local_dofs; r++) {
        d_elementStiffnessMatrix[r].resize(num_local_dofs);
        for(unsigned int c = 0; c < num_local_dofs; c++) {
          d_elementStiffnessMatrix[r][c] = 0;
        }
      }

      d_flowGalWFLinElem->initializeForCurrentElement( d_currElemPtr, d_transportModel );
      d_flowGalWFLinElem->setElementVectors( elementInputVectors );
      d_flowGalWFLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
    }

    void NavierStokesGalWFLinearFEOperator :: postElementOperation() {

      for(unsigned int r = 0; r < d_type0DofIndices.size() ; r++) {
        for(unsigned int dr = 0; dr < 4; dr++) {
          for(unsigned int c = 0; c < d_type0DofIndices.size() ; c++) {
            for(unsigned int dc = 0; dc < 4; dc++) {
              if( dr<3 && dc<3 ){
                d_matrix->addValueByGlobalID( d_type0DofIndices[r][dr], d_type0DofIndices[c][dc], 
                    d_elementStiffnessMatrix[(4*r) + dr][(4*c) + dc] );
              }else{
                d_matrix->addValueByGlobalID( d_type1DofIndices[r], d_type1DofIndices[r], 
                    d_elementStiffnessMatrix[(4*r) + 3][(4*c) + 3] );
              }
            }
          }
        }
      }

      destroyHex27LibMeshElement() ;
    }

    void NavierStokesGalWFLinearFEOperator :: gettype0DofIndicesForCurrentElement(int varId, std::vector<std::vector<size_t> > & dofIds) {
      dofIds.resize(d_currNodes.size());
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds[j]);
      } // end of j
     }

    void NavierStokesGalWFLinearFEOperator :: gettype1DofIndicesForCurrentElement(int varId, std::vector<size_t> & dofIds) {
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_dofMap[varId]->getDOFs(d_currNodes[j].globalID(), dofIds);
      } // end of j
     }

    void NavierStokesGalWFLinearFEOperator :: createHex27LibMeshElement() {
      d_currElemPtr = new ::Hex27;
      for(size_t j = 0; j < d_currNodes.size(); j++) {
        std::vector<double> pt = d_currNodes[j].coord();
        d_currElemPtr->set_node(j) = new ::Node(pt[0], pt[1], pt[2], j);
      }//end for j
    }

    void NavierStokesGalWFLinearFEOperator :: destroyHex27LibMeshElement() {
      for(size_t j = 0; j < d_currElemPtr->n_nodes(); j++) {
        delete (d_currElemPtr->get_node(j));
        d_currElemPtr->set_node(j) = NULL;
      }//end for j
      delete d_currElemPtr;
      d_currElemPtr = NULL;
    }
*/
  }
}//end namespace

