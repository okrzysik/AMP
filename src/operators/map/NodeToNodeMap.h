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

    virtual AMP::LinearAlgebra::Variable::shared_ptr  getInputVariable (int varId = -1) { return d_inpVariable; }
    virtual AMP::LinearAlgebra::Variable::shared_ptr  getOutputVariable () { return d_inpVariable; }

private:

    class Point
    {
       public:
          AMP::Mesh::MeshElementID id;
          double pos[3];
          size_t dof[3];
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

    int dim;
    int DofsPerObj;
    AMP::Discretization::DOFManager::shared_ptr   d_DOFManager1;
    AMP::Discretization::DOFManager::shared_ptr   d_DOFManager2;

    // List of points owned by the current processor for each mesh
    std::vector<Point>                      d_ownedPointsMesh1;
    std::vector<Point>                      d_ownedPointsMesh2;








    // Number of DOFs to send for each mesh for each processor
    std::vector<int>                        d_N_mesh1;
    std::vector<int>                        d_N_mesh2;

    // Displacements for the send for each mesh for each processor
    std::vector<int>                        d_disp_mesh1;
    std::vector<int>                        d_disp_mesh2;

    // List of the local DOFs to send for each mesh for each processor
    std::vector<size_t>                     d_send_DOFs_mesh1;
    std::vector<size_t>                     d_send_DOFs_mesh2;

    // List of the local DOFs to recv for each mesh for each processor
    std::vector<size_t>                     d_recv_DOFs_mesh1;
    std::vector<size_t>                     d_recv_DOFs_mesh2;




    // Buffers for sending/recieving data    

    // Vector of the procs indicating which mesh they own
    // 0: no meshes, 1: mesh1 only, 2: mesh2 only, 3: both meshes
    std::vector<unsigned char>              d_owner;

    // Vector of IDs of elements on the surface
    std::vector<AMP::Mesh::MeshElementID>   d_SurfaceIDs1;
    std::vector<AMP::Mesh::MeshElementID>   d_SurfaceIDs2;

    // Vector of DOFs for the surfaces
    std::vector<size_t>                     d_SurfaceDOFs1;
    std::vector<unsigned int>               d_SurfaceDOFs2;
    
    // Vector to store the send/recv info for the data
    std::vector<double>                     d_SendBuffer;
    std::vector<double>                     d_RecvBuffer;

    // Communicator for the map
    AMP_MPI                                 d_MapComm;

    AMP::LinearAlgebra::Vector::shared_ptr        d_OutputVector;
    AMP::LinearAlgebra::Variable::shared_ptr      d_inpVariable;


    int   d_SendTag;
    int   d_RecvTag;

    void  sendSurface ( boost::shared_ptr<NodeToNodeMapParameters> );
    void  recvSurface ( boost::shared_ptr<NodeToNodeMapParameters> );
    void  sendOrder ( boost::shared_ptr<NodeToNodeMapParameters> );
    void  recvOrder ( boost::shared_ptr<NodeToNodeMapParameters> );
    void  buildSendRecvList ( boost::shared_ptr<NodeToNodeMapParameters> );
    void  finalizeCommunication ( boost::shared_ptr<NodeToNodeMapParameters> );

    class CommInfo
    {
       public:
          AMP::Mesh::MeshElementID   _remId;
          std::list<int>   _procs;
    };


    // Function to create the list of owned points from the iterator over the surface nodes
    std::vector<Point> createOwnedPoints( AMP::Mesh::MeshIterator, AMP::Discretization::DOFManager::shared_ptr );

protected:

};


}
}


#endif
