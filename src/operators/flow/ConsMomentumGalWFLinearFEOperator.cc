#include "ConsMomentumGalWFLinearFEOperator.h"
#include "utils/Utilities.h"

namespace AMP {
  namespace Operator {

    ConsMomentumGalWFLinearFEOperator :: ConsMomentumGalWFLinearFEOperator (
        const boost::shared_ptr<ConsMomentumGalWFLinearFEOperatorParameters> & params)
      : LinearFEOperator (params) {
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

        d_flowGalWFLinElem = boost::dynamic_pointer_cast<ConsMomentumGalWFLinearElement>(d_elemOp);

        AMP_INSIST( ((d_flowGalWFLinElem.get()) != NULL), "d_elemOp is not of type ConsMomentumGalWFLinearElement" );


        d_transportModel = params->d_transportModel ;

        std::string inpVarName = params->d_db->getString("InputVariable");
        d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVarName) );

        std::string outVarName = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVarName) );

        d_inVec.resize(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);
        bool isAttachedToNonlinearOperator = params->d_db->getBoolWithDefault("isAttachedToNonlinearOperator", false);

        if(isAttachedToNonlinearOperator) {
          bool isNonlinearOperatorInitialized = params->d_db->getBoolWithDefault("isNonlinearOperatorInitialized", false);
          if(isNonlinearOperatorInitialized) {
            reset(params);
          } else {
            AMP_ERROR("Not Converted yet");
            //d_matrix = d_MeshAdapter->createMatrix ( d_inpVariable, d_outVariable );
          }
        } else {
          reset(params);
        }
      }

    void ConsMomentumGalWFLinearFEOperator :: preAssembly(const boost::shared_ptr<OperatorParameters>& oparams) 
    {

      boost::shared_ptr<ConsMomentumGalWFLinearFEOperatorParameters> params = boost::dynamic_pointer_cast<ConsMomentumGalWFLinearFEOperatorParameters>(oparams);

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

    void ConsMomentumGalWFLinearFEOperator :: postAssembly()
    {

      d_matrix->makeConsistent ();
    }

    void ConsMomentumGalWFLinearFEOperator :: preElementOperation( const AMP::Mesh::MeshElement & elem ) {
      AMP_ERROR("Not converted yet");
      /*
      unsigned int num_local_dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_dofIndices[i], i);
        num_local_dofs += d_dofIndices[i].size();
      }//end for i

      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[NavierStokes::VELOCITY].resize(num_local_dofs);

      d_numNodesForCurrentElement = elem.numNodes();

      for(unsigned int r = 0; r < d_numNodesForCurrentElement; r++) {
        for(unsigned int d = 0; d < 3; d++) {
          elementInputVectors[NavierStokes::VELOCITY][(3*r) + d] = (d_inVec[NavierStokes::VELOCITY])->getValueByGlobalID( d_dofIndices[d][r] );
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
      */
    }

    void ConsMomentumGalWFLinearFEOperator :: postElementOperation() {

      unsigned int num_local_dofs_per_component = d_dofIndices[0].size();

      for(unsigned int r = 0; r < num_local_dofs_per_component; r++) {
        for(unsigned int dr = 0; dr < 3; dr++) {
          for(unsigned int c = 0; c < num_local_dofs_per_component; c++) {
            for(unsigned int dc = 0; dc < 3; dc++) {
              d_matrix->addValueByGlobalID( d_dofIndices[dr][r], d_dofIndices[dc][c], 
                  d_elementStiffnessMatrix[(3*r) + dr][(3*c) + dc] );
            }
          }
        }
      }

    }

  }
}//end namespace




