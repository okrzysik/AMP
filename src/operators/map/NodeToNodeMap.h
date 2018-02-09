#ifndef included_AMP_NodeToNodeMap
#define included_AMP_NodeToNodeMap


#include <list>
#include <vector>

#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/operators/map/NodeToNodeMapParameters.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"


namespace AMP {
namespace Operator {


/** \brief  A class for mapping a vector from one surface to another where the
 * nodes are aligned.  This routine is bi-directional.  It will exchange data
 * for the aligned nodes.
 */
class NodeToNodeMap : public AMP::Operator::AsyncMapOperator
{
public:
    //! brief  Typedef to identify the parameters class of this operator
    typedef NodeToNodeMapParameters Parameters;

    /** \brief  Returns true if MapType = "NodeToNode"
     * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
     * \return  True iff s == "NodeToNode"
     */
    static bool validMapType( const std::string &s );

    //!  The base tag used in communication.
    static const int CommTagBase = 10000;

    //! Constructor
    explicit NodeToNodeMap( const AMP::shared_ptr<AMP::Operator::OperatorParameters> &params );

    //! Destructor
    virtual ~NodeToNodeMap();

    virtual void applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual void applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                              AMP::LinearAlgebra::Vector::shared_ptr f ) override;

    virtual void setVector( AMP::LinearAlgebra::Vector::shared_ptr p ) override;

    virtual AMP::LinearAlgebra::Vector::shared_ptr getVector();

    // Function to determine if a makeConsistentSet is required
    virtual bool requiresMakeConsistentSet() override;


private:
    class Point
    {
    public:
        AMP::Mesh::MeshElementID id;
        double pos[3];
        int proc;
        Point();
        Point( const Point &rhs );
        bool operator==( const Point &rhs ) const;
        bool operator!=( const Point &rhs ) const;
        bool operator<=( const Point &rhs ) const;
        bool operator>=( const Point &rhs ) const;
        bool operator<( const Point &rhs ) const;
        bool operator>( const Point &rhs ) const;
    };

    // Some basic variables
    AMP::Mesh::MeshIterator d_iterator1;
    AMP::Mesh::MeshIterator d_iterator2;

    // Store the pairs of points that are aligned for each mesh owned by the current processor
    std::vector<std::pair<Point, Point>> d_localPairsMesh1;
    std::vector<std::pair<Point, Point>> d_localPairsMesh2;

    // Variables for communication
    int d_commTag;
    std::vector<int> d_count;
    std::vector<int> d_displ;
    std::vector<double> d_sendBuffer;
    std::vector<double> d_recvBuffer;

    // Other data
    bool d_callMakeConsistentSet;

    // Function to compute the pairs of points for each mesh
    // Note: This function requires global communication across the map comm
    void createPairs( bool requireAllPaired = true );

    // Function to create the list of owned points from the iterator over the surface nodes
    std::vector<Point> createOwnedPoints( const AMP::Mesh::MeshIterator & );

    // Function to create the communication lists
    void buildSendRecvList();

protected:
    int dim;
    int DofsPerObj;
    // Lists of MeshElementIDs to send/recv for each processor (note: orders must match for each
    // processor)
    // We want to construct the lists so that we can do a global or pair-wise communication
    std::vector<AMP::Mesh::MeshElementID> d_sendList;
    std::vector<AMP::Mesh::MeshElementID> d_recvList;
};
} // namespace Operator
} // namespace AMP


#endif
