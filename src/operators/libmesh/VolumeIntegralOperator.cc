#include "operators/libmesh/VolumeIntegralOperator.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "ProfilerApp.h"

#include "vectors/VectorBuilder.h"
#include "matrices/MatrixBuilder.h"

#include <cstring>

namespace AMP {
namespace Operator {


VolumeIntegralOperator::VolumeIntegralOperator(
    const AMP::shared_ptr<VolumeIntegralOperatorParameters> & params) : 
    NonlinearFEOperator(params)
{
    AMP_INSIST( ((params.get()) != NULL), "NULL parameter!" );
    AMP_INSIST( (((params->d_db).get()) != NULL), "NULL database!" );

    d_srcNonlinElem = AMP::dynamic_pointer_cast<SourceNonlinearElement>(d_elemOp);
    AMP_INSIST( ((d_srcNonlinElem.get()) != NULL), "d_elemOp is not of type SourceNonlinearElement" );

    if ( params->d_sourcePhysicsModel ) {
        d_sourcePhysicsModel = params->d_sourcePhysicsModel; 
    } else if (params->d_db->keyExists("SourcePhysicsModel") ) {
        d_sourcePhysicsModel = params->d_sourcePhysicsModel;
    }

    AMP::shared_ptr<AMP::Database> primaryDb = params->d_db->getDatabase("ActiveInputVariables");

    int numPrimaryVariables   =  (params->d_db)->getInteger("Number_Active_Variables");
    int numAuxillaryVariables =  (params->d_db)->getInteger("Number_Auxillary_Variables");

    d_inpVariables.reset(new AMP::LinearAlgebra::MultiVariable("myInpVar"));
    d_auxVariables.reset(new AMP::LinearAlgebra::MultiVariable("myAuxVar"));

    for(int i = 0; i < numPrimaryVariables ; i++) {
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

    for (int var = 0; var < numPrimaryVariables  ; var++)
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


void VolumeIntegralOperator::preAssembly(AMP::LinearAlgebra::Vector::const_shared_ptr u, 
    AMP::LinearAlgebra::Vector::shared_ptr r)
{
    AMP_INSIST( (u != NULL), "NULL Input Vector" );

    AMP::LinearAlgebra::VS_Mesh meshSelector(d_Mesh);
    AMP::LinearAlgebra::Vector::const_shared_ptr meshSubsetPrimary, meshSubsetAuxillary;

    if ( d_inpVariables->numVariables() > 0 )
        meshSubsetPrimary = u->constSelect(meshSelector,d_inpVariables->getName());
    for(size_t var = 0; var < d_inpVariables->numVariables(); var++)
    {
        AMP::LinearAlgebra::Variable::shared_ptr primaryVariable = d_inpVariables->getVariable(var);
        d_inVec[var] = meshSubsetPrimary->constSubsetVectorForVariable( primaryVariable );
        AMP_ASSERT( d_inVec[var] != NULL );
        AMP_ASSERT(d_inVec[var]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);
    }

    if ( d_auxVariables->numVariables() > 0 )
        meshSubsetAuxillary = d_multiAuxPtr->select(meshSelector,d_auxVariables->getName());
    for(size_t var = 0; var < d_auxVariables->numVariables(); var++)
    {
        AMP::LinearAlgebra::Variable::shared_ptr auxillaryVariable = d_auxVariables->getVariable(var);
        d_auxVec[var] = meshSubsetAuxillary->constSubsetVectorForVariable( auxillaryVariable );
        AMP_ASSERT( d_auxVec[var] != NULL );
        AMP_ASSERT(d_auxVec[var]->getUpdateStatus()==AMP::LinearAlgebra::Vector::UNCHANGED);
    }

    // subsetOutputVector is from Operator.h
    d_outVec = this->subsetOutputVector(r);
    d_outVec->zero();

    if(d_inpVariables->numVariables() > 0) {
        d_elementDofMap = d_inVec[0]->getDOFManager();
    } else if(d_auxVariables->numVariables() > 0) {
        d_elementDofMap = d_auxVec[0]->getDOFManager();
    }

    d_nodeDofMap = d_outVec->getDOFManager();

    if ( d_isInputType== "NodalScalar" ) {
        for (unsigned int var = 0; var < d_inpVariables->numVariables(); var++)
        {
            AMP::Mesh::MeshElementID tmp = d_inVec[var]->getDOFManager()->getIterator()->globalID();
            if ( tmp.type() != AMP::Mesh::Vertex )
                AMP_ERROR("Input vector isn't really a NodalScalar");
        }
        for (unsigned int var = 0; var < d_auxVariables->numVariables(); var++)
        {
            AMP::Mesh::MeshElementID tmp = d_auxVec[var]->getDOFManager()->getIterator()->globalID();
            if ( tmp.type() != AMP::Mesh::Vertex )
                AMP_ERROR("aux vector isn't really a NodalScalar");
        } 
    }

}


void VolumeIntegralOperator::preElementOperation(
    const AMP::Mesh::MeshElement & elem )
{
    PROFILE_START("preElementOperation",5);
    d_currNodes = elem.getElements(AMP::Mesh::Vertex);

    std::vector<size_t> elemDofIds;
    d_elementDofMap->getDOFs(elem.globalID(), elemDofIds);

    getNodeDofIndicesForCurrentElement(); 

    std::vector<std::vector<double> > elementInputVectors(d_inpVariables->numVariables());
    std::vector<std::vector<double> > elementAuxVectors(d_auxVariables->numVariables());

    if ( d_isInputType == "IntegrationPointScalar" ) {
        AMP_INSIST(!elemDofIds.empty(),"d_elementDofMap does not contain element, but type is IntegrationPointScalar");
        for (unsigned int var = 0; var < d_inpVariables->numVariables(); var++)
        {
            elementInputVectors[var].resize(elemDofIds.size());
            d_inVec[var]->getValuesByGlobalID( elemDofIds.size(), &elemDofIds[0], &elementInputVectors[var][0] );
        }
        for (unsigned int var = 0; var < d_auxVariables->numVariables(); var++)
        {
            elementAuxVectors[var].resize(elemDofIds.size());
            d_auxVec[var]->getValuesByGlobalID( elemDofIds.size(), &elemDofIds[0], &elementAuxVectors[var][0] );
        }
    } else if( d_isInputType== "NodalScalar" ) {
        AMP_INSIST(elemDofIds.empty(),"d_elementDofMap contains elements, but type is NodalScalar");
        for (unsigned int var = 0; var < d_inpVariables->numVariables(); var++)
        {
            elementInputVectors[var].resize(d_dofIndices.size());
            d_inVec[var]->getValuesByGlobalID( d_dofIndices.size(), &d_dofIndices[0], &elementInputVectors[var][0] );
        }
        for (unsigned int var = 0; var < d_auxVariables->numVariables(); var++)
        {
            elementAuxVectors[var].resize(d_dofIndices.size());
            d_auxVec[var]->getValuesByGlobalID( d_dofIndices.size(), &d_dofIndices[0], &elementAuxVectors[var][0] );
        }
    }

    d_elementOutputVector.resize(d_dofIndices.size());
    // Reinitialize the std::vector to zero (resize does not do this)
    for(unsigned int i = 0; i < d_dofIndices.size(); i++) 
        d_elementOutputVector[i] = 0.0;

    d_srcNonlinElem->initializeForCurrentElement(d_currElemPtrs[d_currElemIdx], d_sourcePhysicsModel);

    d_srcNonlinElem->setElementVectors(elementInputVectors, elementAuxVectors, d_elementOutputVector);
    PROFILE_STOP("preElementOperation",5);
}


void VolumeIntegralOperator::postElementOperation() 
{
    PROFILE_START("postElementOperation",5);
    d_outVec->addValuesByGlobalID( d_dofIndices.size(), &d_dofIndices[0], &d_elementOutputVector[0] );
    PROFILE_STOP("postElementOperation",5);
}


void VolumeIntegralOperator::postAssembly()
{
    d_outVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_ADD );
}


void VolumeIntegralOperator :: init( const AMP::shared_ptr<VolumeIntegralOperatorParameters>& ) 
{
    AMP::Mesh::MeshIterator  el     = d_Mesh->getIterator(AMP::Mesh::Volume, 0);
    d_srcNonlinElem->setElementFlags(d_isInputType);
    for(d_currElemIdx=0; d_currElemIdx<el.size(); ++d_currElemIdx, ++el) {
        d_currNodes = el->getElements(AMP::Mesh::Vertex);
        d_srcNonlinElem->initializeForCurrentElement(d_currElemPtrs[d_currElemIdx], d_sourcePhysicsModel);
    }//end for el
    d_currElemIdx = static_cast<unsigned int>(-1);
}


void VolumeIntegralOperator :: reset( const AMP::shared_ptr<OperatorParameters>& )
{
    d_outVec.reset();
}


AMP::shared_ptr<OperatorParameters> VolumeIntegralOperator::getJacobianParameters(
    const AMP::shared_ptr<AMP::LinearAlgebra::Vector>& u ) 
{
    (void) u;
    AMP::shared_ptr<AMP::InputDatabase> tmp_db( new AMP::InputDatabase("Dummy"));
    tmp_db->putString("name", "VolumeIntegralOperator");
    AMP::shared_ptr<VolumeIntegralOperatorParameters> outParams( new VolumeIntegralOperatorParameters(tmp_db));

    outParams->d_sourcePhysicsModel= d_sourcePhysicsModel;
    outParams->d_pVector = u;
    return outParams;
}


/*
AMP::shared_ptr<AMP::LinearAlgebra::Matrix> VolumeIntegralOperator::getLinearizedVolumeIntegralOperator(
   const AMP::shared_ptr<OperatorParameters>& params)            
{    
   AMP::shared_ptr<VolumeIntegralOperatorParameters> inParams = AMP::dynamic_pointer_cast<VolumeIntegralOperatorParameters>(params);
   const AMP::shared_ptr<AMP::LinearAlgebra::Vector> u = inParams->d_pVector;

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
    std::vector<size_t> dofs;
    for(unsigned int j = 0; j < d_currNodes.size(); j++) {
        d_nodeDofMap->getDOFs(d_currNodes[j].globalID(),dofs);
        AMP_ASSERT(dofs.size()==1);
        d_dofIndices[j] = dofs[0];
    }// end of j
}


}
}//end namespace



