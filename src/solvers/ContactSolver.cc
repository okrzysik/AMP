
#include "ContactSolver.h"
#include "operators/ContactResidualCorrection.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
  namespace Solver {

    void ContactSolver :: solve(boost::shared_ptr<AMP::LinearAlgebra::Vector>, boost::shared_ptr<AMP::LinearAlgebra::Vector>  u) {
      boost::shared_ptr<AMP::Operator::ContactResidualCorrection> op = boost::dynamic_pointer_cast<
        AMP::Operator::ContactResidualCorrection>(d_pOperator);

      AMP::LinearAlgebra::Variable::shared_ptr masterVariable = op->getMasterVariable();
      AMP::LinearAlgebra::Variable::shared_ptr slaveVariable = op->getSlaveVariable();

      std::vector<AMP::Mesh::MeshElementID> masterNodes = op->getMasterNodes();
      std::vector<AMP::Mesh::MeshElementID> slaveNodes = op->getSlaveNodes();

      std::vector<std::vector<unsigned int> > dofs = op->getDofs();

      AMP::LinearAlgebra::Vector::shared_ptr uMaster = u->subsetVectorForVariable(masterVariable);
      AMP::LinearAlgebra::Vector::shared_ptr uSlave = u->subsetVectorForVariable(slaveVariable);

      AMP::Discretization::DOFManager::shared_ptr master_dof_map = uMaster->getDOFManager();
      AMP::Discretization::DOFManager::shared_ptr slave_dof_map = uSlave->getDOFManager();

      for(size_t i = 0; i < masterNodes.size(); i++) {
        std::vector<size_t> masterGlobalIds;
        std::vector<size_t> slaveGlobalIds;
        master_dof_map->getDOFs(masterNodes[i], masterGlobalIds);
        slave_dof_map->getDOFs(slaveNodes[i], slaveGlobalIds);
        for(size_t j = 0; j < dofs[i].size(); j++) {
          double masterVal = uMaster->getLocalValueByGlobalID( masterGlobalIds[dofs[i][j]] );
          uSlave->setLocalValueByGlobalID(slaveGlobalIds[dofs[i][j]], masterVal);
        }//end for j
      }//end for i
    }

  }
}



