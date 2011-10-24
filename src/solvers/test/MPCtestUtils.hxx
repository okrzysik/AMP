
void getMasterAndSlaveNodes(std::vector<unsigned int> & masterNodes, std::vector<unsigned int> & slaveNodes,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter) {
  AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd;
  AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd;

  bnd = masterMeshAdapter->beginOwnedBoundary( 2 );
  end_bnd = masterMeshAdapter->endOwnedBoundary( 2 );

  unsigned int numSurfaceNodes = (end_bnd.getIterator() - bnd.getIterator());

  masterNodes.resize(numSurfaceNodes);
  slaveNodes.resize(numSurfaceNodes);

  for(int i = 0; bnd != end_bnd; i++, bnd++) {
    masterNodes[i] = bnd->globalID();
  }

  bnd = slaveMeshAdapter->beginOwnedBoundary( 1 );
  end_bnd = slaveMeshAdapter->endOwnedBoundary( 1 );

  for(int i = 0; bnd != end_bnd; i++, bnd++) {
    slaveNodes[i] = bnd->globalID();
  }

  for(unsigned int i = 0; i < numSurfaceNodes; i++) {
    std::cout<<"master["<<i<<"] = "<<masterNodes[i]<<std::endl;
    std::cout<<"slave["<<i<<"] = "<<slaveNodes[i]<<std::endl;
  }

  std::cout<<"# Interface Nodes = "<<numSurfaceNodes<<std::endl;
}

void setSlaveToZero(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Vector::shared_ptr vec,
    AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter, std::vector<unsigned int> slaveNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd =  slaveMeshAdapter->getNode( slaveNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);
    }//end for j
  }//end for i
}

void addSlaveToMaster(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    std::vector<unsigned int> slaveNodes, std::vector<unsigned int> masterNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd = slaveMeshAdapter->getNode( slaveNodes[i] );
    AMP::Mesh::LibMeshNode masterNd = masterMeshAdapter->getNode( masterNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    std::vector<unsigned int> masterGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      double slaveVal = slaveVec->getLocalValueByGlobalID(slaveGlobalIds[j]);
      masterVec->addLocalValueByGlobalID(masterGlobalIds[j], slaveVal);
    }//end for j
  }//end for i
}

void copyMasterToSlave(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
    AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
    std::vector<unsigned int> slaveNodes, std::vector<unsigned int> masterNodes) {

  AMP::LinearAlgebra::Vector::shared_ptr slaveVec = vec->subsetVectorForVariable(slaveVar);
  AMP::LinearAlgebra::Vector::shared_ptr masterVec = vec->subsetVectorForVariable(masterVar);

  AMP::Mesh::DOFMap::shared_ptr slave_dof_map = slaveMeshAdapter->getDOFMap(slaveVar);
  AMP::Mesh::DOFMap::shared_ptr master_dof_map = masterMeshAdapter->getDOFMap(masterVar);

  std::vector<unsigned int> dofs(3);
  dofs[0] = 0; dofs[1] = 1; dofs[2] = 2;

  for(size_t i = 0; i < slaveNodes.size(); i++) {
    AMP::Mesh::LibMeshNode slaveNd = slaveMeshAdapter->getNode( slaveNodes[i] );
    AMP::Mesh::LibMeshNode masterNd = masterMeshAdapter->getNode( masterNodes[i] );
    std::vector<unsigned int> slaveGlobalIds;
    std::vector<unsigned int> masterGlobalIds;
    slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, dofs);
    master_dof_map->getDOFs(masterNd, masterGlobalIds, dofs);
    for(size_t j = 0; j < dofs.size(); j++) {
      double masterVal = masterVec->getLocalValueByGlobalID(masterGlobalIds[j]);
      slaveVec->setLocalValueByGlobalID(slaveGlobalIds[j], masterVal);
    }//end for j
  }//end for i
}




