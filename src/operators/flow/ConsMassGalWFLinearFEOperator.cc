#include "ConsMassGalWFLinearFEOperator.h"
#include "utils/Utilities.h"


namespace AMP {
namespace Operator {

ConsMassGalWFLinearFEOperator :: ConsMassGalWFLinearFEOperator (
    const boost::shared_ptr<ConsMassGalWFLinearFEOperatorParameters> & params)
    : LinearFEOperator (params) 
{
        AMP_INSIST( ((params.get()) != NULL), "NULL parameter" );

        d_flowGalWFLinElem = boost::dynamic_pointer_cast<ConsMassGalWFLinearElement>(d_elemOp);

        AMP_INSIST( ((d_flowGalWFLinElem.get()) != NULL), "d_elemOp is not of type ConsMassGalWFLinearElement" );


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
            AMP_ERROR("Not converted, createMatrix has been moved to MatrixBuilder and required the right and left vectors");
            //d_matrix = d_MeshAdapter->createMatrix ( d_inpVariable, d_outVariable );
          }
        } else {
          reset(params);
        }
}


void ConsMassGalWFLinearFEOperator :: preAssembly(const boost::shared_ptr<OperatorParameters>& oparams) 
{

      boost::shared_ptr<ConsMassGalWFLinearFEOperatorParameters> params = boost::dynamic_pointer_cast<ConsMassGalWFLinearFEOperatorParameters>(oparams);

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


void ConsMassGalWFLinearFEOperator :: postAssembly()
{
      d_matrix->makeConsistent ();
}


void ConsMassGalWFLinearFEOperator :: preElementOperation( const AMP::Mesh::MeshElement & elem ) 
{
AMP_ERROR("Not converted yet"); /*
      unsigned int num_local_u_dofs = 0;
      for(unsigned int i = 0; i < 3; i++) {
        (dof_maps[0])->getDOFs (elem, d_dofIndices0[i], i);
        num_local_u_dofs += d_dofIndices0[i].size();
      }//end for i

      //    unsigned int dofIds=0;
      (dof_maps[1])->getDOFs (elem, d_dofIndices1);
      unsigned int num_local_p_dofs = d_dofIndices1.size();

      std::vector<std::vector<double> > elementInputVectors(NavierStokes::TOTAL_NUMBER_OF_VARIABLES);

      elementInputVectors[NavierStokes::PRESSURE].resize(num_local_p_dofs);

      d_numNodesForCurrentElement = elem.numNodes();

      //   for(unsigned int r = 0; r < num_local_p_dofs; r++) {
      //   elementInputVectors[NavierStokes::PRESSURE][r] = (d_inVec[NavierStokes::PRESSURE])->getValueByGlobalID( d_dofIndices1[r] );
      //   }

      d_elementStiffnessMatrix.resize(num_local_u_dofs);
      for(unsigned int r = 0; r < num_local_u_dofs; r++) {
        d_elementStiffnessMatrix[r].resize(num_local_p_dofs);
        for(unsigned int c = 0; c < num_local_p_dofs; c++) {
          d_elementStiffnessMatrix[r][c] = 0;
        }
      }

      const ::Elem* elemPtr = &(elem.getElem());

      d_flowGalWFLinElem->initializeForCurrentElement( elemPtr, d_transportModel );
      d_flowGalWFLinElem->setElementVectors( elementInputVectors );
      d_flowGalWFLinElem->setElementStiffnessMatrix( d_elementStiffnessMatrix );
  */
}


void ConsMassGalWFLinearFEOperator :: postElementOperation() {

    unsigned int num_local_u_dofs = d_dofIndices0[0].size();
    unsigned int num_local_p_dofs = d_dofIndices1.size();

    for(unsigned int r = 0; r < num_local_u_dofs; r++) {
        for(unsigned int dr = 0; dr < 3; dr++) {
            for(unsigned int c = 0; c < num_local_p_dofs; c++) {
                d_matrix->addValueByGlobalID( d_dofIndices0[dr][r], d_dofIndices1[c], 
                    d_elementStiffnessMatrix[(3*r) + dr][c] );
            }
        }
    }
}


}//end namespace
}



