
#include "VolumeIntegralOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

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

        d_numPrimaryVariables   =  (params->d_db)->getInteger("Number_Active_Variables");
        d_numAuxillaryVariables =  (params->d_db)->getInteger("Number_Auxillary_Variables");

        d_isInputType = (params->d_db)->getStringWithDefault("InputVariableType", "NodalScalar");

        d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));

        for(unsigned int i = 0; i < d_numPrimaryVariables ; i++) {
            AMP::LinearAlgebra::Variable::shared_ptr dummyVar;
            d_inpVariables->add(dummyVar);
        }
                
        d_inVec.resize(d_numPrimaryVariables);

        if(d_numAuxillaryVariables>0){
            AMP_INSIST( (((params->d_auxVec).get()) != NULL), "NULL Auxillary Vector!" );}
                d_multiAuxPtr = params->d_auxVec;
        d_auxVec.resize(d_numAuxillaryVariables);

                d_activeVariableNames.resize(d_numPrimaryVariables);

        for (unsigned int var = 0; var < d_numPrimaryVariables  ; var++)
        {
            char key[100];
            if(d_isInputType == "IntegrationPointScalar"){
                sprintf(key, "ActiveVariable_%d", (int)var);
                std::string varName = primaryDb->getString(key);
                                d_activeVariableNames[var] = varName;
                AMP::LinearAlgebra::Variable::shared_ptr inpVar(new AMP::LinearAlgebra::Variable(varName));
                d_inpVariables->setVariable(var, inpVar);
            }else if(d_isInputType== "NodalScalar"){
                sprintf(key, "ActiveVariable_%d", (int)var);
                std::string varName = primaryDb->getString(key);
                                d_activeVariableNames[var] = varName;
                AMP::LinearAlgebra::Variable::shared_ptr inpVar(new AMP::LinearAlgebra::Variable(varName));
                d_inpVariables->setVariable(var, inpVar);
            }
        }

//            unsigned int nPts = (params->d_db)->getIntegerWithDefault("number_integration_points",8);
//            d_inpVariable.reset(new AMP::Mesh::RunTimeIntegrationPointVariable(inpVar, nPts));

        std::string outVar = params->d_db->getString("OutputVariable");
        d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar)); 

        d_bMatrixAndVectorsCloned=false;

        init(params);
    }

    void VolumeIntegralOperator::preAssembly(
            const boost::shared_ptr<AMP::LinearAlgebra::Vector> &u, boost::shared_ptr<AMP::LinearAlgebra::Vector> &r)
    {
        AMP_INSIST( (u != NULL), "NULL Input Vector" );

        for (unsigned int var = 0; var < d_numPrimaryVariables ; var++)
        {
            AMP::LinearAlgebra::Variable::shared_ptr primaryVariable = d_inpVariables->getVariable(var);
        

            d_inVec[var] = u->subsetVectorForVariable( primaryVariable );
            AMP_ASSERT( d_inVec[var] != NULL );
            (d_inVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }
        for (unsigned int var = 0; var < d_numAuxillaryVariables ; var++)
        {
            AMP::LinearAlgebra::Variable::shared_ptr auxillaryVariable = d_auxVariables->getVariable(var);

            d_auxVec[var] = d_multiAuxPtr->subsetVectorForVariable( auxillaryVariable );
            (d_auxVec[var])->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        }

        d_outVec = r->subsetVectorForVariable(d_outVariable);
        d_outVec->zero();

    }

    void VolumeIntegralOperator::postAssembly()
    {
        d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
    }

    void VolumeIntegralOperator::preElementOperation(
            const AMP::Mesh::MeshElement & elem )
    {
    AMP_ERROR("Not converted yet");
    /*
        // unsigned int num_local_Dofs = 0;
        dof_maps[0]->getDOFs (elem, d_inpDofIndices);
        // num_local_Dofs = d_inpDofIndices.size();

        AMP::Mesh::DOFMap::shared_ptr temp_dof_maps;
        temp_dof_maps = d_MeshAdapter->getDOFMap( d_outVariable );
        temp_dof_maps->getDOFs (elem, d_outDofIndices);

        std::vector<std::vector<double> > elementInputVectors(d_numPrimaryVariables);
        std::vector<std::vector<double> > elementAuxVectors(d_numAuxillaryVariables);

        d_numNodesForCurrentElement = elem.numNodes();

        if(d_isInputType == "IntegrationPointScalar"){
            d_numDofsForCurrentElement =  dof_maps[0]->numDOFsPerObject(); 
            for (unsigned int var = 0; var < d_numPrimaryVariables; var++)
            {
                elementInputVectors[var].resize(d_numDofsForCurrentElement);
                for (unsigned int i = 0; i < d_numDofsForCurrentElement  ; i++) {
                    elementInputVectors[var][i] =  d_inVec[var]->getValueByGlobalID( dof_maps[0]->getGlobalID(elem.globalID(), i) );
                }
            }
            for (unsigned int var = 0; var < d_numAuxillaryVariables; var++)
            {
                elementAuxVectors[var].resize(d_numDofsForCurrentElement);
                for (unsigned int i = 0; i < d_numDofsForCurrentElement  ; i++) {
                    elementAuxVectors[var][i] =  d_auxVec[var]->getValueByGlobalID( dof_maps[0]->getGlobalID(elem.globalID(), i) );
                }
            }
        }else if(d_isInputType== "NodalScalar"){
            for (unsigned int var = 0; var < d_numPrimaryVariables; var++)
            {
                elementInputVectors[var].resize(d_numNodesForCurrentElement);
                for (unsigned int i = 0; i < d_numNodesForCurrentElement  ; i++) {
                    elementInputVectors[var][i] =  d_inVec[var]->getValueByGlobalID( d_inpDofIndices[i] );
                }
            }
            for (unsigned int var = 0; var < d_numAuxillaryVariables; var++)
            {
                elementAuxVectors[var].resize(d_numNodesForCurrentElement);
                for (unsigned int i = 0; i < d_numNodesForCurrentElement  ; i++) {
                    elementAuxVectors[var][i] =  d_auxVec[var]->getValueByGlobalID( d_inpDofIndices[i] );
                }
            }
        }

        d_elementOutputVector.resize(d_numNodesForCurrentElement);
        for(unsigned int i = 0; i < d_numNodesForCurrentElement; i++) {
            d_elementOutputVector[i] = 0.0;
        }

        const ::Elem* elemPtr = &(elem.getElem());

        d_srcNonlinElem->initializeForCurrentElement(elemPtr,d_sourcePhysicsModel);

        d_srcNonlinElem->setElementVectors(elementInputVectors, elementAuxVectors, d_elementOutputVector);
    */
    }

    void VolumeIntegralOperator::postElementOperation() 
    {

        for (unsigned int i = 0; i < d_numNodesForCurrentElement ; i++) {
            d_outVec->addValueByGlobalID(d_outDofIndices[i], d_elementOutputVector[i]);
        }
    }

    void VolumeIntegralOperator :: init(const boost::shared_ptr<VolumeIntegralOperatorParameters>& params) 
    {
        int ghostWidth = 0;
        AMP::Mesh::MeshIterator  el     = d_Mesh->getIterator(AMP::Mesh::Volume, ghostWidth);
        AMP::Mesh::MeshIterator  end_el = el.end();

        d_srcNonlinElem->setElementFlags(d_isInputType);
        for( ; el != end_el; ++el) {
/*
          const ::Elem* elemPtr = &(el->getElem());
          d_srcNonlinElem->initializeForCurrentElement( elemPtr, d_sourcePhysicsModel);
*/
        }//end for el
        
    }

    void VolumeIntegralOperator :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      AMP_ERROR("Not converted yet");
      /*
         AMP::Mesh::DOFMap::shared_ptr dof_maps;

      //    dof_maps = d_MeshAdapter->getDOFMap( d_inpVariables );

      // AMP::Mesh::MeshManager::Adapter::ElementIterator  el = d_MeshAdapter->beginElement();
      // AMP::Mesh::MeshManager::Adapter::ElementIterator  end_el = d_MeshAdapter->endElement();

      boost::shared_ptr<VolumeIntegralOperatorParameters> myParams =
      boost::dynamic_pointer_cast<VolumeIntegralOperatorParameters>(params); 

      AMP_INSIST( ((myParams.get()) != NULL), "Null parameter!" );

      d_outVec.reset();
      */
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

    boost::shared_ptr<AMP::LinearAlgebra::Matrix> VolumeIntegralOperator::getLinearizedVolumeIntegralOperator(
        const boost::shared_ptr<OperatorParameters>& params)            
    {    
      AMP_ERROR("Not converted yet");
      /*    
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

            */      
            return d_pDiagonalMatrix;  
    }


}
}//end namespace

