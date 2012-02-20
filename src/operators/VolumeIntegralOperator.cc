
#include "VolumeIntegralOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

#include <cstring>

namespace AMP {
  namespace Operator {
    VolumeIntegralOperator::VolumeIntegralOperator(
        const boost::shared_ptr<VolumeIntegralOperatorParameters> & params) : NonlinearFEOperator(params)
    {
      AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
      AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

      d_srcNonlinElem = boost::dynamic_pointer_cast<SourceNonlinearElement>(d_elemOp);
      AMP_INSIST( ((d_srcNonlinElem.get()) != NULL), "d_elemOp is not of type SourceNonlinearElement" );

      if ( params->d_sourcePhysicsModel ) {
        d_sourcePhysicsModel = params->d_sourcePhysicsModel; 
      } else if (params->d_db->keyExists("SourcePhysicsModel") ) {
        d_sourcePhysicsModel = params->d_sourcePhysicsModel;
      }

      boost::shared_ptr<AMP::Database> primaryDb = params->d_db->getDatabase("ActiveInputVariables");

      int numPrimaryVariables   =  (params->d_db)->getInteger("Number_Active_Variables");
      int numAuxillaryVariables =  (params->d_db)->getInteger("Number_Auxillary_Variables");

      d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
      d_auxVariables.reset(new AMP::LinearAlgebra::MultiVariable("myAuxVar"));

      for(unsigned int i = 0; i < numPrimaryVariables ; i++) {
        AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
        d_inpVariables->add(dummyVar);
      }

      d_inVec.resize(numPrimaryVariables);

      AMP_INSIST( ( numAuxillaryVariables == 0), "Verify this works before using it; the Interface to SourcePhysicsModel.h does not appear to be complete." );
      if(numAuxillaryVariables>0){
        AMP_INSIST( (((params->d_auxVec).get()) != NULL), "NULL Auxillary Vector!" );
      }
      d_multiAuxPtr = params->d_auxVec;
      d_auxVec.resize(numAuxillaryVariables);

      for (unsigned int var = 0; var < numPrimaryVariables  ; var++)
      {
        char key[100];
        sprintf(key, "ActiveVariable_%d", (int)var);
        std::string varName = primaryDb->getString(key);
        AMP::LinearAlgebra::Variable::shared_ptr inpVar(new AMP::LinearAlgebra::Variable(varName));
        d_inpVariables->setVariable(var, inpVar);
      }

      std::string outVar = params->d_db->getString("OutputVariable");
      d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar)); 

      d_isInputType = params->d_db->getStringWithDefault("InputVariableType", "IntegrationPointScalar");

      //d_bMatrixAndVectorsCloned=false;

      init(params);
    }

    void VolumeIntegralOperator::preAssembly(const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, 
        boost::shared_ptr<AMP::LinearAlgebra::Vector> &r)
    {
      AMP_INSIST( (u != NULL), "NULL Input Vector" );

      for(size_t var = 0; var < d_inpVariables->numVariables(); var++)
      {
        AMP::LinearAlgebra::Variable::shared_ptr primaryVariable = d_inpVariables->getVariable(var);
        d_inVec[var] = u->subsetVectorForVariable( primaryVariable );
        AMP_ASSERT( d_inVec[var] != NULL );
        (d_inVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }

      for(size_t var = 0; var < d_auxVariables->numVariables(); var++)
      {
        AMP::LinearAlgebra::Variable::shared_ptr auxillaryVariable = d_auxVariables->getVariable(var);
        d_auxVec[var] = d_multiAuxPtr->subsetVectorForVariable( auxillaryVariable );
        (d_auxVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
      }

      d_outVec = r->subsetVectorForVariable(d_outVariable);
      d_outVec->zero();

      if(d_inpVariables->numVariables() > 0) {
        d_elementDofMap = d_inVec[0]->getDOFManager();
      } else if(d_auxVariables->numVariables() > 0) {
        d_elementDofMap = d_auxVec[0]->getDOFManager();
      }

      d_nodeDofMap = d_outVec->getDOFManager();
    }

    void VolumeIntegralOperator::preElementOperation(
        const AMP::Mesh::MeshElement & elem )
    {
      d_currNodes = elem.getElements(AMP::Mesh::Vertex);

      createCurrentLibMeshElement();

      std::vector<size_t> elemDofIds;
      d_elementDofMap->getDOFs(elem.globalID(), elemDofIds);

      getNodeDofIndicesForCurrentElement(); 

      std::vector<std::vector<double> > elementInputVectors(d_inpVariables->numVariables());
      std::vector<std::vector<double> > elementAuxVectors(d_auxVariables->numVariables());

      if(d_isInputType == "IntegrationPointScalar"){
        for (unsigned int var = 0; var < d_inpVariables->numVariables(); var++)
        {
          elementInputVectors[var].resize(elemDofIds.size());
          for(size_t i = 0; i < elemDofIds.size(); i++) {
            elementInputVectors[var][i] =  d_inVec[var]->getValueByGlobalID( elemDofIds[i] );
          }
        }
        for (unsigned int var = 0; var < d_auxVariables->numVariables(); var++)
        {
          elementAuxVectors[var].resize(elemDofIds.size());
          for (size_t i = 0; i < elemDofIds.size(); i++) {
            elementAuxVectors[var][i] =  d_auxVec[var]->getValueByGlobalID( elemDofIds[i] );
          }
        }
      }else if(d_isInputType== "NodalScalar"){
        for (unsigned int var = 0; var < d_inpVariables->numVariables(); var++)
        {
          elementInputVectors[var].resize(d_dofIndices.size());
          for (size_t i = 0; i < d_dofIndices.size(); i++) {
            elementInputVectors[var][i] =  d_inVec[var]->getValueByGlobalID( d_dofIndices[i][0] );
          }
        }
        for (unsigned int var = 0; var < d_auxVariables->numVariables(); var++)
        {
          elementAuxVectors[var].resize(d_dofIndices.size());
          for(size_t i = 0; i < d_dofIndices.size(); i++) {
            elementAuxVectors[var][i] =  d_auxVec[var]->getValueByGlobalID( d_dofIndices[i][0] );
          }
        } 
      }

      d_elementOutputVector.resize(d_dofIndices.size(), 0.0);

      d_srcNonlinElem->initializeForCurrentElement(d_currElemPtr,d_sourcePhysicsModel);

      d_srcNonlinElem->setElementVectors(elementInputVectors, elementAuxVectors, d_elementOutputVector);
    }

    void VolumeIntegralOperator::postElementOperation() 
    {
      for (size_t i = 0; i < d_dofIndices.size(); i++) {
        d_outVec->addValueByGlobalID(d_dofIndices[i][0], d_elementOutputVector[i]);
      }
      destroyCurrentLibMeshElement();
    }

    void VolumeIntegralOperator::postAssembly()
    {
      d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    void VolumeIntegralOperator :: init(const boost::shared_ptr<VolumeIntegralOperatorParameters>& params) 
    {
      AMP::Mesh::MeshIterator  el     = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
      AMP::Mesh::MeshIterator  end_el = el.end();
      d_srcNonlinElem->setElementFlags(d_isInputType);
      for( ; el != end_el; ++el) {
        d_currNodes = el->getElements(AMP::Mesh::Vertex);
        createCurrentLibMeshElement();
        d_srcNonlinElem->initializeForCurrentElement(d_currElemPtr, d_sourcePhysicsModel);
        destroyCurrentLibMeshElement();
      }//end for el
    }

    void VolumeIntegralOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      d_outVec.reset();
    }

    boost::shared_ptr<OperatorParameters> VolumeIntegralOperator::getJacobianParameters(
        const boost::shared_ptr<AMP::LinearAlgebra::Vector>& u) 
    {
      (void) u;
      boost::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase("Dummy"));
      tmp_db->putString("name", "VolumeIntegralOperator");
      boost::shared_ptr<VolumeIntegralOperatorParameters> outParams( new VolumeIntegralOperatorParameters(tmp_db));

      outParams->d_sourcePhysicsModel= d_sourcePhysicsModel;
      outParams->d_pVector = u;
      return outParams;
    }

    /*
       boost::shared_ptr<AMP::LinearAlgebra::Matrix> VolumeIntegralOperator::getLinearizedVolumeIntegralOperator(
       const boost::shared_ptr<OperatorParameters>& params)            
       {    
       boost::shared_ptr<VolumeIntegralOperatorParameters> inParams = boost::dynamic_pointer_cast<VolumeIntegralOperatorParameters>(params);
       const boost::shared_ptr<AMP::LinearAlgebra::Vector> u = inParams->d_pVector;

       if (!d_bMatrixAndVectorsCloned)
       {
       d_pDiagonalMatrix = d_MeshAdapter->createMatrix(this->getInputVariable(0),this->getOutputVariable());
       d_pDiagonalVector = d_MeshAdapter->createVector(this->getOutputVariable());
       d_pNullVector     = d_MeshAdapter->createVector(this->getOutputVariable());

       AMP::pout << "in the loop" << std::endl;
       d_bMatrixAndVectorsCloned=true;
       }

       AMP::pout << "d_pDiagonalVector.get() = " << d_pDiagonalVector.get() << std::endl;
       this->apply(d_pNullVector,u,d_pDiagonalVector,1.0,0.0);
       d_pDiagonalMatrix->setDiagonal(d_pDiagonalVector);

       return d_pDiagonalMatrix;  
       }
       */

    void VolumeIntegralOperator::getNodeDofIndicesForCurrentElement() {
      d_dofIndices.resize(d_currNodes.size());
      for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_nodeDofMap->getDOFs(d_currNodes[j].globalID(), d_dofIndices[j]);
      }// end of j
    }

  }
}//end namespace



