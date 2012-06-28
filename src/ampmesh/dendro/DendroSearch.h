
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
    bool verbose;
    AMP::AMP_MPI globalComm;
    AMP::Mesh::Mesh::shared_ptr meshAdapter;
    int rank, npes;
    double minCoords[3];
    double maxCoords[3];
    double ScalingFactor[3];
    std::vector<ot::TreeNode> nodeList;
    std::vector<int> stIdxList;
    std::vector<int> rankList;
    std::vector<int> elemIdList;
    std::vector<AMP::Mesh::MeshElement> localElemArr;
    std::vector<ot::TreeNode> mins;
    unsigned int BoxLevel;

    int numLocalPts;
    std::vector<double> foundPts;
    std::vector<int> sendCnts, sendDisps, recvCnts, recvDisps;

    void setupDendro();
};

#endif // DENDRO_SEARCH
