#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"


void getMasterAndSlaveNodes(
    std::vector<AMP::Mesh::MeshElementID> & masterNodes, 
    std::vector<AMP::Mesh::MeshElementID> & slaveNodes,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter) 
{

  AMP::Mesh::MeshIterator bnd;
  AMP::Mesh::MeshIterator end_bnd;

  bnd = masterMeshAdapter->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
  end_bnd = bnd.end();

  unsigned int numSurfaceNodes = bnd.size();

  masterNodes.resize(numSurfaceNodes);
  slaveNodes.resize(numSurfaceNodes);

  for(int i = 0; bnd != end_bnd; i++, bnd++) {
    masterNodes[i] = bnd->globalID();
  }

  bnd = slaveMeshAdapter->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
  end_bnd = bnd.end();

  for(int i = 0; bnd != end_bnd; i++, bnd++) {
    slaveNodes[i] = bnd->globalID();
  }

  //for(unsigned int i = 0; i < numSurfaceNodes; i++) {
  //  std::cout<<"master["<<i<<"] = "<<masterNodes[i]<<std::endl;
  //  std::cout<<"slave["<<i<<"] = "<<slaveNodes[i]<<std::endl;
  //}

  std::cout<<"# Interface Nodes = "<<numSurfaceNodes<<std::endl;
}

void setSlaveToZero(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Vector::shared_ptr vec,
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter, std::vector<AMP::Mesh::MeshElementID> slaveNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::Discretization::DOFManager::shared_ptr slave_dof_map = slaveVec->getDOFManager();

  std::vector<size_t> slaveGlobalIds;
  for(size_t i = 0; i < slaveNodes.size(); i++) {
    slave_dof_map->getDOFs( slaveNodes[i], slaveGlobalIds);
    for(size_t j = 0; j < slaveGlobalIds.size(); j++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);
    }//end for j
  }//end for i
}

void addSlaveToMaster(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    std::vector<AMP::Mesh::MeshElementID> slaveNodes, std::vector<AMP::Mesh::MeshElementID> masterNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  AMP::Discretization::DOFManager::shared_ptr  slave_dof_map = slaveVec->getDOFManager();
  AMP::Discretization::DOFManager::shared_ptr  master_dof_map = masterVec->getDOFManager();

  std::vector<size_t> slaveGlobalIds, masterGlobalIds;
  for(size_t i = 0; i < slaveNodes.size(); i++) {
    slave_dof_map->getDOFs( slaveNodes[i], slaveGlobalIds );
    master_dof_map->getDOFs( masterNodes[i], masterGlobalIds );
    for(size_t j = 0; j < slaveGlobalIds.size(); j++) {
      double slaveVal = slaveVec->getLocalValueByGlobalID(slaveGlobalIds[j]);
      masterVec->addLocalValueByGlobalID(masterGlobalIds[j], slaveVal);
    }//end for j
  }//end for i
}

void copyMasterToSlave(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    std::vector<AMP::Mesh::MeshElementID> slaveNodes, std::vector<AMP::Mesh::MeshElementID> masterNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  AMP::Discretization::DOFManager::shared_ptr  slave_dof_map = slaveVec->getDOFManager();
  AMP::Discretization::DOFManager::shared_ptr  master_dof_map = masterVec->getDOFManager();

  std::vector<size_t> slaveGlobalIds, masterGlobalIds;
  for(size_t i = 0; i < slaveNodes.size(); i++) {
    slave_dof_map->getDOFs( slaveNodes[i], slaveGlobalIds );
    master_dof_map->getDOFs( masterNodes[i], masterGlobalIds );
    for(size_t j = 0; j < slaveGlobalIds.size(); j++) {
      double masterVal = masterVec->getLocalValueByGlobalID(masterGlobalIds[j]);
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], masterVal);
    }//end for j
  }//end for i
}




