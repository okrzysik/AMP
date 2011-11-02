#if 0
//#ifdef USE_AMP_MESH
#include "WriteSolutionToFile.h"

void printSolution(AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter,
    AMP::LinearAlgebra::Vector::shared_ptr solVec, std::string exeName) {

  AMP::Mesh::DOFMap::shared_ptr dof_map = meshAdapter->getDOFMap(solVec->getVariable());

  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  nd = meshAdapter->beginOwnedNode();
  AMP::Mesh::MeshManager::Adapter::OwnedNodeIterator  end_nd = meshAdapter->endOwnedNode();

  std::string fname = "results_" + exeName + ".txt";
  FILE* fp = fopen(fname.c_str(), "w");

  fprintf(fp, "%s\n\n", exeName.c_str());

  fprintf(fp, "x, y, z,   u,  v,  w\n\n");

  for(; nd != end_nd; ++nd ) {
    std::vector<unsigned int> ndGlobalIds;
    std::vector<unsigned int> dofIds;
    dof_map->getDOFs(*nd, ndGlobalIds, dofIds);

    fprintf(fp, "%lf, %lf, %lf,    ", (nd->x()), (nd->y()), (nd->z()));

    for(int i = 0; i < 3; i++) {
      double val = solVec->getLocalValueByGlobalID(ndGlobalIds[i]);
      fprintf(fp, "%.13lf, ", val);
    }//end for i

    fprintf(fp, " \n");
  }//end for nd

  fclose(fp);
}

#endif

