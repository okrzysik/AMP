
#ifndef DENDRO_SEARCH
#define DENDRO_SEARCH

#include "utils/AMPManager.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/hex8_element_t.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

#include "sys/sys.h"
#include "par/parUtils.h"
#include "binOps/binUtils.h"
#include "oct/TreeNode.h"
#include "oct/octUtils.h"
#include "oct/nodeAndValues.h"
#include "oct/nodeAndRanks.h"
#include "dendro.h"


class DendroSearch {
  public:
    DendroSearch(AMP::AMP_MPI comm, AMP::Mesh::Mesh::shared_ptr mesh);

    void interpolate(AMP::LinearAlgebra::Vector::shared_ptr vectorField, const unsigned int dofsPerNode,
        const std::vector<double> & pts, std::vector<double> & results, std::vector<bool> & foundPt);

    void search(const std::vector<double> & pts);

    void interpolate(AMP::LinearAlgebra::Vector::shared_ptr vectorField, const unsigned int dofsPerNode,
        std::vector<double> & results, std::vector<bool> & foundPt);

  private:
    bool d_verbose;
    AMP::AMP_MPI d_globalComm;
    AMP::Mesh::Mesh::shared_ptr d_meshAdapter;
    int d_rank, d_npes;
    double d_minCoords[3];
    double d_maxCoords[3];
    double d_scalingFactor[3];
    std::vector<ot::TreeNode> d_nodeList;
    std::vector<int> d_stIdxList;
    std::vector<int> d_rankList;
    std::vector<int> d_elemIdList;
    std::vector<AMP::Mesh::MeshElement> d_localElemArr;
    std::vector<ot::TreeNode> d_mins;
    unsigned int d_boxLevel;
    int d_numLocalPts;
    std::vector<double> d_foundPts;
    std::vector<int> d_sendCnts, d_sendDisps, d_recvCnts, d_recvDisps;

    void setupDendro();
    void setupDSforSearch();
    void createLocalMeshElementArray();
};

#endif // DENDRO_SEARCH


