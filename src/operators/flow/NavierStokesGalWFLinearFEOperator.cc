
#include "NavierStokesGalWFLinearFEOperator.h"
#include "utils/Utilities.h"

#if 0
//This file has not been converted!

namespace AMP {
  namespace Operator {

    NavierStokesGalWFLinearFEOperator :: NavierStokesGalWFLinearFEOperator (
        const boost::shared_ptr<NavierStokesGalWFLinearFEOperatorParameters> & params)
      : LinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

        d_flowGalWFLinElem = boost::dynamic_pointer_cast<NavierStokesGalWFLinearElement>(d_elemOp);

        AMP_INSIST( ((d_flowGalWFLinElem.get()) != NULL), "d_elemOp is not of type NavierStokesGalWFLinearElement" );


        d_transportModel = params->d_transportModel ;

        d_inpVariable.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
        d_outVariable.reset(new AMP::LinearAlgebra::MultiVariable("myOutVar"));
        for(unsigned int i = 0; i < NavierStokes::TOTAL_NUMBER_OF_VARIABLES; i++) {
          AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
          d_inpVariable->add(dummyVar);
          d_outVariable->add(dummyVar);
        }//end for i

        AMP_INSIST( params->d_db->keyExists("ActiveInputVariables"), "key not found" );
        boost::shared_ptr<AMP::Database> activeInpVar_db = params->d_db->getDatabase("ActiveInputVariables");

        AMP_INSIST(activeInpVar_db->keyExists("VELOCITY"), "VELOCITY must be active");
        AMP_INSIST(activeInpVar_db->keyExists("PRESSURE"), "VELOCITY must be active");

        std::string tempVarName = activeInpVar_db->getString("VELOCITY");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 3>(tempVarName, d_MeshAdapter) ); 
        d_inpVariable->setVariable(NavierStokes::VELOCITY, tempVar);
        d_outVariable->setVariable(NavierStokes::VELOCITY, tempVar);

        tempVarName = activeInpVar_db->getString("PRESSURE");
        AMP::LinearAlgebra::Variable::shared_ptr tempVar2(new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable, 1>(tempVarName, d_MeshAdapter) ); 
        d_inpVariable->setVariable(NavierStokes::PRESSURE, tempVar2);
        d_outVariable->setVariable(NavierStokes::PRESSURE, tempVar2);

        bool isAttachedToNonlinearOperator = params->d_db->getBoolWithDefault("isAttachedToNonlinearOperator", false);

        if(isAttachedToNonlinearOperator) {
          bool isNonlinearOperatorInitialized = params->d_db->getBoolWithDefault("isNonlinearOperatorInitialized", false);
          if(isNonlinearOperatorInitialized) {
            reset(params);
          } else {
            d_matrix = d_MeshAdapter->createMatrix ( d_inpVariable, d_outVariable );
          }
        } else {
          reset(params);
        }
      }

    void NavierStokesGalWFLinearFEOperator :: preAssembly(const boost::shared_ptr<OperatorParameters>& oparams) 
    {

      boost::shared_ptr<NavierStokesGalWFLinearFEOperatorParameters> params = boost::dynamic_pointer_cast<NavierStokesGalWFLinearFEOperatorParameters>(oparams);

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

    void NavierStokesGalWFLinearFEOperator :: preElementOperation( const AMP::Mesh::MeshManager::Adapter::Element & elem ) {
      unsigned int num_local_dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_dofIndices[i], i);
        num_local_dofs += d_dofIndices[i].size();
      }//end for i

      unsigned int num_local_dofs1 = 0;
      (dof_maps[1])->getDOFs (elem, d_dofIndices1);
      num_local_dofs1 += d_dofIndices1.size();

      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[NavierStokes::VELOCITY].resize(num_local_dofs);
      elementInputVectors[NavierStokes::PRESSURE].resize(num_local_dofs1);
      if(d_inVec[NavierStokes::TEMPERATURE].get() != NULL) {
        elementInputVectors[NavierStokes::TEMPERATURE].resize(num_local_dofs1);
      }

      d_numNodesForCurrentElement = elem.numNodes();

      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->getValueByGlobalID( d_dofIndices[d][r] );
        }

        elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->getValueByGlobalID( d_dofIndices1[r] );

        if(d_inVec[NavierStokes::TEMPERATURE].get() != NULL) {
          elementInputVectors[NavierStokes::TEMPERATURE][r] = (d_inVec[NavierStokes::TEMPERATURE])->getValueByGlobalID( d_dofIndices1[r] );
        }
      }


      d_elementStiffnessMatrix.resize(num_local_dofs);
      for(unsigned int r = 0; r < num_local_dofs; r++) {
        d_elementStiffnessMatrix[r].resize(num_local_dofs);
        for(unsigned int c = 0; c < num_local_dofs; c++) {
          d_elementStiffnessMatrix[r][c] = 0;
        }
      }

      const ::Elem* elemPtr = &(elem.getElem());

      d_flowGalWFLinElem->initializeForCurrentElement( elemPtr, d_transportModel );
      d_flowGalWFLinElem->setElementVectors( elementInputVectors );
      d_flowGalWFLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
    }

    void NavierStokesGalWFLinearFEOperator :: postElementOperation() {

      unsigned int num_local_dofs_per_component = d_dofIndices[0].size();

      for(unsigned int r = 0; r < num_local_dofs_per_component; r++) {
        for(unsigned int dr = 0; dr < 4; dr++) {
          for(unsigned int c = 0; c < num_local_dofs_per_component; c++) {
            for(unsigned int dc = 0; dc < 4; dc++) {
              if( dr<3 && dc<3 ){
                d_matrix->addValueByGlobalID( d_dofIndices[dr][r], d_dofIndices[dc][c], 
                    d_elementStiffnessMatrix[(4*r) + dr][(4*c) + dc] );
              }else{
                d_matrix->addValueByGlobalID( d_dofIndices1[r], d_dofIndices1[c], 
                    d_elementStiffnessMatrix[(4*r) + dr][(4*c) + dc] );
              }
            }
          }
        }
      }

    }

  }
}//end namespace

#endif


