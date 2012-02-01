#ifndef included_AMP_NodeToNodeMap
#define included_AMP_NodeToNodeMap


#include <vector>
#include <list>

#include "operators/map/AsyncMapOperator.h"
#include "operators/map/NodeToNodeMapParameters.h"

#include "vectors/Vector.h"
#include "ampmesh/Mesh.h"
#include "utils/AMP_MPI.h"
#include "discretization/DOF_Manager.h"


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
    typedef  NodeToNodeMapParameters   Parameters;

    /** \brief  Returns true if MapType = "NodeToNode"
      * \param[in] s  A string extracted from the MapType line in a MeshToMeshMap db
      * \return  True iff s == "NodeToNode"
      */
    static bool  validMapType ( const std::string &s );

    //!  The base tag used in communication.
    enum { CommTagBase = 10000 };

    //! Constructor
      NodeToNodeMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> & params );

    //! Destructor
    virtual ~NodeToNodeMap ();

    virtual void applyStart(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

    virtual void applyFinish(const AMP::LinearAlgebra::Vector::shared_ptr &f, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr  &r, const double a = -1.0, const double b = 1.0);

    virtual void  setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p );

    virtual AMP::LinearAlgebra::Variable::shared_ptr  getInputVariable () { return d_inpVariable; }
    virtual AMP::LinearAlgebra::Variable::shared_ptr  getOutputVariable () { return d_inpVariable; }

private:

    class Point
    {
       public:
          AMP::Mesh::MeshElementID id;
          double pos[3];
          size_t dof[8];
          int proc;
          Point ();
          Point ( const Point &rhs );
          bool operator == ( const Point &rhs ) const;
          bool operator != ( const Point &rhs ) const;
          bool operator <= ( const Point &rhs ) const;
          bool operator >= ( const Point &rhs ) const;
          bool operator <  ( const Point &rhs ) const;
          bool operator >  ( const Point &rhs ) const;
    };

    // Some basic variables
    int dim;
    int DofsPerObj;
    AMP::Mesh::MeshIterator  d_iterator1;
    AMP::Mesh::MeshIterator  d_iterator2;
    AMP::Discretization::DOFManager::shared_ptr   d_DOFManager;

    // Store the pairs of points that are aligned for each mesh owned by the current processor
    std::vector< std::pair<Point,Point> >   d_localPairsMesh1;
    std::vector< std::pair<Point,Point> >   d_localPairsMesh2;

    // Lists of DOFs to send/recv for each processor (note: orders must match for each processor)
    // We want to construct the lists so that we can do a global or pair-wise communication
    std::vector<size_t>                     d_sendList;
    std::vector<size_t>                     d_recvList;

    // Variables for communication
    int                                     d_commTag;
    std::vector<int>                        d_count;
    std::vector<int>                        d_displ;
    std::vector<double>                     d_sendBuffer;
    std::vector<double>                     d_recvBuffer;

    // Other data
    //AMP::LinearAlgebra::Vector::shared_ptr        d_OutputVector;
    AMP::LinearAlgebra::Variable::shared_ptr      d_inpVariable;
    bool d_callMakeConsistentSet;

    // Function to compute the pairs of points for each mesh
    // Note: This function requires global communication across the map comm
    void createPairs( bool requireAllPaired=true );

    // Function to create the list of owned points from the iterator over the surface nodes
    std::vector<Point> createOwnedPoints( AMP::Mesh::MeshIterator, AMP::Discretization::DOFManager::shared_ptr );

    // Function to create the communication lists
    void  buildSendRecvList( );

    // Function to determine if a makeConsistentSet is required
    virtual bool requiresMakeConsistentSet();

protected:

};


}
}


#endif
