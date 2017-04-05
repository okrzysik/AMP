
#ifndef DENDRO_SEARCH
#define DENDRO_SEARCH

#include "utils/AMP_MPI.h"
#include "utils/Utilities.h"

#include "vectors/Vector.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/hex8_element_t.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "binOps/binUtils.h"
#include "dendro.h"
#include "oct/TreeNode.h"
#include "oct/nodeAndRanks.h"
#include "oct/nodeAndValues.h"
#include "oct/octUtils.h"
#include "par/parUtils.h"


namespace AMP {
namespace Mesh {


class DendroSearch
{
public:
    struct ProjectOnBoundaryData;

    enum SearchStatus { NotFound = 0, Found, FoundNotOnBoundary, FoundOnBoundary };

    enum TimingType {
        Setup = 0,
        CoarseSearch,
        FineSearch,
        Interpolation,
        ProjectionOnBoundaryID,
        numTimingTypes
    };

    DendroSearch( AMP::Mesh::Mesh::shared_ptr mesh,
                  bool verbose,
                  std::ostream &oStream = std::cout );

    ~DendroSearch()
    {
        for ( size_t i = 0; i < d_volume_elements.size(); ++i ) {
            delete ( d_volume_elements[i] );
        } // end i
        d_volume_elements.clear();
    }

    void searchAndInterpolate( AMP::AMP_MPI comm,
                               AMP::LinearAlgebra::Vector::const_shared_ptr vectorField,
                               const unsigned int dofsPerNode,
                               const std::vector<double> &pts,
                               std::vector<double> &results,
                               std::vector<bool> &foundPt );

    void search( AMP::AMP_MPI comm, const std::vector<double> &pts );

    void interpolate( AMP::AMP_MPI comm,
                      AMP::LinearAlgebra::Vector::const_shared_ptr vectorField,
                      const unsigned int dofsPerNode,
                      std::vector<double> &results,
                      std::vector<bool> &foundPt );

    // deprecated
    void projectOnBoundaryID( AMP::AMP_MPI comm,
                              const int boundaryID,
                              std::vector<AMP::Mesh::MeshElementID> &faceVerticesGlobalIDs,
                              std::vector<double> &shiftGlobalCoords,
                              std::vector<double> &projectionLocalCoordsOnGeomType::Face,
                              std::vector<int> &flags );

    void projectOnBoundaryID( AMP::AMP_MPI comm,
                              const int boundaryID,
                              std::vector<AMP::Mesh::MeshElementID> &faceVerticesGlobalIDs,
                              std::vector<double> &shiftGlobalCoords,
                              std::vector<double> &projectionLocalCoordsOnGeomType::Face,
                              std::vector<int> &flags,
                              std::vector<AMP::Mesh::MeshElementID> &volumeGlobalIDs,
                              std::vector<size_t> &faceLocalIndices );

    void setTolerance( double tolerance );

    void reportTiming( size_t n, TimingType const *timingTypes, double *timingMeasurements );

private:
    AMP::Mesh::Mesh::shared_ptr d_meshAdapter;
    std::vector<AMP::Mesh::MeshElement> d_localElems;
    std::vector<hex8_element_t *> d_volume_elements;
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
    std::ostream &d_oStream;
    std::vector<double> d_timingMeasurements;
    double d_tolerance;

    void setupDSforSearch();
};

struct DendroSearch::ProjectOnBoundaryData {
    size_t d_PointLocalID;
    SearchStatus d_SearchStatus;
    AMP::Mesh::MeshElementID d_GeomType::FaceVerticesIDs[4];
    double d_ProjectionLocalCoordsOnGeomType::Face[2];
    double d_ShiftGlobalCoords[3];
    size_t d_GeomType::FaceLocalIndex;
    AMP::Mesh::MeshElementID d_GeomType::VolumeID;
};
}
}

#endif // DENDRO_SEARCH
