
#include "ContactResidualCorrection.h"

namespace AMP {
  namespace Operator {

    void ContactResidualCorrection :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double , const double ) {
      AMP::LinearAlgebra::Vector::shared_ptr rMaster = r->subsetVectorForVariable(d_masterVariable);
      AMP::LinearAlgebra::Vector::shared_ptr rSlave = r->subsetVectorForVariable(d_slaveVariable);

      AMP::Mesh::DOFMap::shared_ptr master_dof_map = d_MeshAdapter->getDOFMap(d_masterVariable);
      AMP::Mesh::DOFMap::shared_ptr slave_dof_map = d_slaveMeshAdapter->getDOFMap(d_slaveVariable);

      for(size_t i = 0; i < d_masterNodes.size(); i++) {
        AMP::Mesh::LibMeshNode masterNd =  d_MeshAdapter->getNode( d_masterNodes[i] );
        AMP::Mesh::LibMeshNode slaveNd =  d_slaveMeshAdapter->getNode( d_slaveNodes[i] );
        std::vector<unsigned int> masterGlobalIds;
        std::vector<unsigned int> slaveGlobalIds;
        master_dof_map->getDOFs(masterNd, masterGlobalIds, d_dofs[i]);
        slave_dof_map->getDOFs(slaveNd, slaveGlobalIds, d_dofs[i]);
        for(size_t j = 0; j < d_dofs[i].size(); j++) {
          double slaveVal = rSlave->getLocalValueByGlobalID( slaveGlobalIds[j] );
          double masterVal = rMaster->getLocalValueByGlobalID( masterGlobalIds[j] );

          rMaster->addLocalValueByGlobalID(masterGlobalIds[j], slaveVal);

          masterVal = rMaster->getLocalValueByGlobalID( masterGlobalIds[j] );

          rSlave->setLocalValueByGlobalID(slaveGlobalIds[j], 0.0);

          slaveVal = rSlave->getLocalValueByGlobalID( slaveGlobalIds[j] );
        }//end for j
      }//end for i
    }

  }
}


