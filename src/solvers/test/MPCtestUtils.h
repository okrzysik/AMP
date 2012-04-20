
#ifndef included_MPCtestUtils
#define included_MPCtestUtils


void getMasterAndSlaveNodes(std::vector<unsigned int> & masterNodes, std::vector<unsigned int> & slaveNodes,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter);

void setSlaveToZero(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Vector::shared_ptr vec,
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter, std::vector<unsigned int> slaveNodes);

void addSlaveToMaster(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    std::vector<unsigned int> slaveNodes, std::vector<unsigned int> masterNodes);

void copyMasterToSlave(AMP::LinearAlgebra::Variable::shared_ptr slaveVar, AMP::LinearAlgebra::Variable::shared_ptr masterVar, 
    AMP::LinearAlgebra::Vector::shared_ptr vec, AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter,
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter,
    std::vector<unsigned int> slaveNodes, std::vector<unsigned int> masterNodes);

#include "MPCtestUtils.hxx"

#endif

