
#ifndef DENDRO_SEARCH
#define DENDRO_SEARCH

#include "utils/AMPManager.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/PIO.h"

#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/hex8_element_t.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>

#include "par/parUtils.h"
#include "binOps/binUtils.h"
#include "oct/TreeNode.h"
#include "oct/octUtils.h"
#include "oct/nodeAndValues.h"
#include "oct/nodeAndRanks.h"
#include "dendro.h"


namespace AMP {
namespace Mesh {


class DendroSearch {
  public:
    struct ProjectOnBoundaryData; 

    enum SearchStatus { NotFound = 0, Found, FoundNotOnBoundary, FoundOnBoundary };

    enum TimingType { Setup = 0, CoarseSearch, FineSearch, Interpolation, ProjectionOnBoundaryID, numTimingTypes };

    DendroSearch(AMP::Mesh::Mesh::shared_ptr mesh, bool verbose = true, std::ostream & oStream = std::cout);

    void searchAndInterpolate(AMP::AMP_MPI comm, AMP::LinearAlgebra::Vector::const_shared_ptr vectorField, const unsigned int dofsPerNode,
        const std::vector<double> & pts, std::vector<double> & results, std::vector<bool> & foundPt);

    void search(AMP::AMP_MPI comm, const std::vector<double> & pts);

    void interpolate(AMP::AMP_MPI comm, AMP::LinearAlgebra::Vector::const_shared_ptr vectorField, const unsigned int dofsPerNode,
        std::vector<double> & results, std::vector<bool> & foundPt);

    void projectOnBoundaryID(AMP::AMP_MPI comm, const int boundaryID, std::vector<AMP::Mesh::MeshElementID> & faceVerticesGlobalIDs, 
        std::vector<double> & shiftGlobalCoords, std::vector<double> & projectionLocalCoordsOnFace, std::vector<int> & flags);

    void setTolerance(double tolerance);

    void reportTiming(size_t n, TimingType const * timingTypes, double * timingMeasurements);

  private:
    AMP::Mesh::Mesh::shared_ptr d_meshAdapter;
    std::vector<AMP::Mesh::MeshElement> d_localElemArr;
    std::vector<ot::TreeNode> d_nodeList;
    std::vector<ot::TreeNode> d_mins;
    std::vector<double> d_minCoords;
    std::vector<double> d_scalingFactor;
    std::vector<double> d_foundPts;
    std::vector<int> d_stIdxList;
    std::vector<int> d_rankList;
    std::vector<int> d_elemIdList;
    std::vector<int> d_sendCnts;
    std::vector<int> d_sendDisps;
    std::vector<int> d_recvCnts;
    std::vector<int> d_recvDisps;
    unsigned int d_boxLevel;
    int d_numLocalPts;

    bool d_verbose;
    std::ostream & d_oStream;
    std::vector<double> d_timingMeasurements;
  
    std::vector<hex8_element_t> d_volume_elements;
    double d_tolerance;

    void setupDSforSearch();
    void createLocalMeshElementArray();
};

struct DendroSearch::ProjectOnBoundaryData {
  size_t d_PointLocalID;
  SearchStatus d_SearchStatus;
  AMP::Mesh::MeshElementID d_FaceVerticesIDs[4];
  double d_ProjectionLocalCoordsOnFace[2]; 
  double d_ShiftGlobalCoords[3];
};


}
}

#endif // DENDRO_SEARCH


