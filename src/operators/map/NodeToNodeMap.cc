#include "operators/map/NodeToNodeMap.h"
#include "operators/map/NodeToNodeMapParameters.h"
#include "ampmesh/MeshElement.h"
#include "discretization/NodalVariable.h"
#include "discretization/simpleDOF_Manager.h"

#include <set>


namespace AMP {
namespace Operator {


template <class T>
static T* getPtr( std::vector<T> &x ) {
    if ( x.size()== 0 )
        return NULL;
    return &x[0];
}



/********************************************************
* Constructor                                           *
********************************************************/
NodeToNodeMap::NodeToNodeMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> & params )
       : AMP::Operator::AsyncMapOperator ( params )
{
    // Cast the params appropriately
    d_OutputVector = AMP::LinearAlgebra::Vector::shared_ptr ();
    AMP_ASSERT ( params );
    NodeToNodeMapParameters &Params =
                        *(boost::dynamic_pointer_cast<NodeToNodeMapParameters> ( params ) );

    // Set class members
    d_MapComm = Params.d_MapComm;
    dim = -1;
    if ( d_mesh1.get() != NULL )
        dim = d_mesh1->getDim();
    dim = d_MapComm.maxReduce(dim);
    AMP_INSIST(dim<=3,"Node to Node map only works up to 3d (see Point)");
    DofsPerObj = Params.d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj<=3,"Node to Node map only works for <= 3 DOFs per node (see Point)");
    d_DOFManager = Params.d_DOFManager;
    d_commTag = Params.d_commTag;


    // Create a nodal variable and DOFManager (this should be moved out of here)
    std::string variableName = Params.d_db->getString("VariableName");
    d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::Discretization::NodalVariable(DofsPerObj,variableName) );

    // Create the element iterators
    if ( d_mesh1.get() != NULL )
        d_iterator1 = d_mesh1->getIDsetIterator(AMP::Mesh::Vertex,Params.d_BoundaryID1,0);
    if ( d_mesh2.get() != NULL )
        d_iterator2 = d_mesh2->getIDsetIterator(AMP::Mesh::Vertex,Params.d_BoundaryID2,0);

    // Create the pairs of points that are aligned
    createPairs();

    // Now we need to create the DOF lists and setup the communication
    buildSendRecvList();

    // Determine the number of pair-wise communications
    size_t numPartners = 0;
    for (size_t i=0; i<d_count.size(); i++) {
        if ( d_count[i] > 0 )
            numPartners++;
    }
    reserveRequests( 2*numPartners );

    // Finished setting up the communication
}


/********************************************************
* De-constructor                                        *
********************************************************/
NodeToNodeMap::~NodeToNodeMap ()
{
}


/********************************************************
* Check if the string matches a NodeToNode map          *
********************************************************/
bool NodeToNodeMap::validMapType ( const std::string &t )
{
    if ( t == "NodeToNode" )
        return true;
    return false;
}


/********************************************************
* Set the vector                                        *
********************************************************/
void  NodeToNodeMap::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p )
{
    d_OutputVector = p->subsetVectorForVariable ( d_inpVariable );
    AMP_INSIST ( d_OutputVector , "setVector received bogus stuff" );
}


/********************************************************
* Start the communication                               *
********************************************************/
void NodeToNodeMap::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
                              const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                                    AMP::LinearAlgebra::Vector::shared_ptr & ,
                              const double ,
                              const double )
{
    // Subset the vector for the variable
    AMP::LinearAlgebra::Vector::shared_ptr   curPhysics = u->subsetVectorForVariable ( d_inpVariable );
    AMP_INSIST( curPhysics , "apply received bogus stuff" );
    AMP_INSIST( curPhysics->getVariable()->DOFsPerObject()==(size_t)DofsPerObj, "vector has different # of DOFs per node than map" );
    AMP_INSIST( curPhysics->getDOFManager()==d_DOFManager,"The DOF Manager that created the vector must match the one for the map" );

    // Get the DOFs to send
    std::vector<int> DOFs(d_sendList.size(),0);
    for (size_t i=0; i<d_sendList.size(); i++)
        DOFs[i] = (int) d_sendList[i];
    curPhysics->getValuesByGlobalID( d_sendList.size(), getPtr( DOFs ), getPtr( d_sendBuffer ) );

    // Start the communication
    std::vector<MPI_Request>::iterator  curReq = beginRequests();
    for (int i=0; i<d_MapComm.getSize(); i++) {
        if ( i==d_MapComm.getRank() ) {
            // Perform a local copy
            for (int j=d_displ[i]; j<d_displ[i]+d_count[i]; j++)
                d_recvBuffer[j] = d_sendBuffer[j];
        } else if ( d_count[i] > 0 ) {
            // Start asyncronous communication
            *curReq = d_MapComm.Isend( &d_sendBuffer[d_displ[i]], d_count[i], i, d_commTag );
            curReq++;
            *curReq = d_MapComm.Irecv( &d_recvBuffer[d_displ[i]], d_count[i], i, d_commTag );
            curReq++;
        }
    }
}


/********************************************************
* Finish the communication and fill the vector          *
********************************************************/
void NodeToNodeMap::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
                              const AMP::LinearAlgebra::Vector::shared_ptr &   ,
                                    AMP::LinearAlgebra::Vector::shared_ptr &r  ,
                              const double ,
                              const double )
{
    // Wait to recieve all data
    waitForAllRequests();
    // bool  copyToOutput = false;

    // Get the vector to store the DOFs
    AMP::LinearAlgebra::Vector::shared_ptr  curPhysics = d_OutputVector;
    //    if ( r )
    //    {
    //      curPhysics = r->subsetVectorForVariable ( d_inpVariable );
    //    }
    //    else if ( d_OutputVector )
    //    {
    //      copyToOutput = false;
    //      curPhysics = d_OutputVector;
    //    }
    AMP_INSIST( curPhysics->getDOFManager()==d_DOFManager,"The DOF Manager that created the vector must match the one for the map" );

    // Store the DOFs
    std::vector<int> DOFs(d_recvList.size(),0);
    for (size_t i=0; i<d_recvList.size(); i++)
        DOFs[i] = (int) d_recvList[i];
    curPhysics->setValuesByGlobalID( d_recvList.size(),  getPtr( DOFs ), getPtr( d_recvBuffer ) );

    // Update ghost cells
    curPhysics->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
}


/****************************************************************
* Function to compute the send/recv lists                       *
* We need to group the DOFs to send/recieve by processor and    *
* ensure that both processors will agree on the order.          *
****************************************************************/
void NodeToNodeMap::buildSendRecvList( )
{
    // First, count the total number of DOFs to send/recv for each processor
    int commSize = d_MapComm.getSize();
    d_count = std::vector<int>(commSize,0);
    d_displ = std::vector<int>(commSize,0);
    for (size_t i=0; i<d_localPairsMesh1.size(); i++) {
        int rank = d_localPairsMesh1[i].second.proc;
        d_count[rank]+=DofsPerObj;
    }
    for (size_t i=0; i<d_localPairsMesh2.size(); i++) {
        int rank = d_localPairsMesh2[i].second.proc;
        d_count[rank]+=DofsPerObj;
    }
    d_displ = std::vector<int>(commSize,0);
    for (int i=1; i<commSize; i++)
        d_displ[i] = d_displ[i-1] + d_count[i-1];
    int N_tot = d_displ[commSize-1] + d_count[commSize-1];
    // Create the send/recv lists and remote DOF lists for each processor
    std::vector< std::vector<size_t> > send_DOFs(commSize);
    std::vector< std::vector<size_t> > recv_DOFs(commSize);
    std::vector< std::vector<size_t> > remote_DOFs(commSize);
    for (int i=0; i<commSize; i++) {
        send_DOFs[i].reserve(d_count[i]);
        recv_DOFs[i].reserve(d_count[i]);
        remote_DOFs[i].reserve(d_count[i]);
    }
    for (size_t i=0; i<d_localPairsMesh1.size(); i++) {
        int rank = d_localPairsMesh1[i].second.proc;
        for (int j=0; j<DofsPerObj; j++) {
            send_DOFs[rank].push_back(d_localPairsMesh1[i].first.dof[j]);
            recv_DOFs[rank].push_back(d_localPairsMesh1[i].first.dof[j]);
            remote_DOFs[rank].push_back(d_localPairsMesh1[i].second.dof[j]);
        }
    }
    for (size_t i=0; i<d_localPairsMesh2.size(); i++) {
        int rank = d_localPairsMesh2[i].second.proc;
        for (int j=0; j<DofsPerObj; j++) {
            send_DOFs[rank].push_back(d_localPairsMesh2[i].first.dof[j]);
            recv_DOFs[rank].push_back(d_localPairsMesh2[i].first.dof[j]);
            remote_DOFs[rank].push_back(d_localPairsMesh2[i].second.dof[j]);
        }
    }
    // Sort the send/recv lists by the sending processor's DOF
    for (int i=0; i<commSize; i++) {
        AMP::Utilities::quicksort(send_DOFs[i]);
        AMP::Utilities::quicksort(remote_DOFs[i],recv_DOFs[i]);
    }
    // Create the final lists and allocate space for the send/recv
    d_sendList = std::vector<size_t>(N_tot,-1);
    d_recvList = std::vector<size_t>(N_tot,-1);
    d_sendBuffer = std::vector<double>(N_tot,0);
    d_recvBuffer = std::vector<double>(N_tot,0);
    for (int i=0; i<commSize; i++) {
        int k = d_displ[i];
        for (size_t j=0; j<send_DOFs[i].size(); j++) {
            d_sendList[k+j] = send_DOFs[i][j];
            d_recvList[k+j] = recv_DOFs[i][j];
        }
    }
}


/****************************************************************
* Function to compute the pairs of points for each mesh         *
* We have the option of requiring if all points are matched.    *
****************************************************************/
void NodeToNodeMap::createPairs( bool requireAllPaired )
{
    d_localPairsMesh1.resize(0);
    d_localPairsMesh2.resize(0);

    // For each mesh, get the list of points owned by the current processor
    std::vector<Point>  ownedPointsMesh1;
    std::vector<Point>  ownedPointsMesh2;
    if ( d_mesh1.get() != NULL )
        ownedPointsMesh1 = createOwnedPoints( d_iterator1, d_DOFManager );
    if ( d_mesh2.get() != NULL )
        ownedPointsMesh2 = createOwnedPoints( d_iterator2, d_DOFManager );


    // Send the list of points on mesh1 to all processors
    int commSize = d_MapComm.getSize();
    int send_cnt = (int) ownedPointsMesh1.size();
    std::vector<int> recv_cnt(commSize,0);
    std::vector<int> recv_disp(commSize,0);
    d_MapComm.allGather(send_cnt,&recv_cnt[0]);
    for (int i=1; i<commSize; i++)
        recv_disp[i] = recv_disp[i-1]+recv_cnt[i-1];
    int N_recv_tot = recv_disp[commSize-1] + recv_cnt[commSize-1];
    std::vector<Point> surfacePts = std::vector<Point>(N_recv_tot);
    d_MapComm.allGather( getPtr(ownedPointsMesh1), send_cnt, &surfacePts[0], &recv_cnt[0], &recv_disp[0], true );

    // Sort the points for fast searching
    AMP::Utilities::quicksort(surfacePts);

    // Find the points in mesh1 that align with the points owned by the current processor on mesh2
    d_localPairsMesh2.reserve(ownedPointsMesh2.size());
    for (size_t i=0; i<ownedPointsMesh2.size(); i++) {
        // Search for the point
        int index = AMP::Utilities::findfirst(surfacePts,ownedPointsMesh2[i]);
        // Check if the point was found
        bool found = index>=0 && index<(int)surfacePts.size();
        if ( found ) { if ( surfacePts[index] != ownedPointsMesh2[i] ) { found = false; } }
        // Add the pair to the list
        if ( found ) {
            std::pair<Point,Point> pair(ownedPointsMesh2[i],surfacePts[index]);
            d_localPairsMesh2.push_back(pair);
        } else if ( requireAllPaired ) {
            AMP_ERROR("All points are required to be paired, and some points were not found");
        }
    }
    surfacePts.clear();


    // Send the list of points on mesh2 to all processors
    send_cnt = (int) ownedPointsMesh2.size();
    d_MapComm.allGather(send_cnt,&recv_cnt[0]);
    for (int i=1; i<commSize; i++)
        recv_disp[i] = recv_disp[i-1]+recv_cnt[i-1];
    N_recv_tot = recv_disp[commSize-1] + recv_cnt[commSize-1];
    surfacePts = std::vector<Point>(N_recv_tot);
    d_MapComm.allGather( getPtr(ownedPointsMesh2), send_cnt, &surfacePts[0], &recv_cnt[0], &recv_disp[0], true );

    // Sort the points for fast searching
    AMP::Utilities::quicksort(surfacePts);

    // Find the points in mesh1 that align with the points owned by the current processor on mesh2
    d_localPairsMesh1.reserve(ownedPointsMesh1.size());
    for (size_t i=0; i<ownedPointsMesh1.size(); i++) {
        // Search for the point
        int index = AMP::Utilities::findfirst(surfacePts,ownedPointsMesh1[i]);
        // Check if the point was found
        bool found = index>=0 && index<(int)surfacePts.size();
        if ( found ) { if ( surfacePts[index] != ownedPointsMesh1[i] ) { found = false; } }
        // Add the pair to the list
        if ( found ) {
            std::pair<Point,Point> pair(ownedPointsMesh1[i],surfacePts[index]);
            d_localPairsMesh1.push_back(pair);
        } else if ( requireAllPaired ) {
            AMP_ERROR("All points are required to be paired, and some points were not found");
        }
    }
    surfacePts.clear();

    // Syncronize and return
    d_MapComm.barrier();
}


/********************************************************
* Function to create the list of owned points from the  *
* iterator over the surface nodes                       *
********************************************************/
std::vector<NodeToNodeMap::Point> NodeToNodeMap::createOwnedPoints( 
    AMP::Mesh::MeshIterator iterator, AMP::Discretization::DOFManager::shared_ptr DOFManager )
{
    // Create the list of points for each node
    std::vector<Point> surfacePts(iterator.size());
    AMP::Mesh::MeshIterator cur = iterator.begin();
    int rank = d_MapComm.getRank();
    std::vector<unsigned int> dofs(DofsPerObj,-1);
    for (size_t i=0; i<surfacePts.size(); i++) {
        // Get the properties of the current element
        AMP::Mesh::MeshElementID id = cur->globalID();
        std::vector<double> pos = cur->coord();
        if ( !id.is_local() )
            continue;
        // Create the point
        Point temp;
        temp.id = cur->globalID();
        DOFManager->getDOFs(temp.id,dofs);
        AMP_INSIST((int)dofs.size()==DofsPerObj,"The specified number of DOFs per object does not match the DOFManager");
        for (int j=0; j<dim; j++)
            temp.pos[j] = pos[j];
        for (int j=0; j<DofsPerObj; j++)
            temp.dof[j] = dofs[j];
        temp.proc = rank;
        surfacePts[i] = temp;
        ++cur;
    }
    // Sort the points
    AMP::Utilities::quicksort(surfacePts);
    return surfacePts;
}


/********************************************************
* Constructors for Point                                *
********************************************************/
NodeToNodeMap::Point::Point ()
{
    id = AMP::Mesh::MeshElementID();
    proc = -1;
    for ( size_t i=0; i!=3; i++)
      pos[i] = 0.0;
    for ( size_t i=0; i!=3; i++)
      dof[i] = -1;
}
NodeToNodeMap::Point::Point ( const Point &rhs )
{
    id = rhs.id;
    proc = rhs.proc;
    for ( size_t i=0; i!=3; i++)
      pos[i] = rhs.pos[i];
    for ( size_t i=0; i!=3; i++)
      dof[i] = rhs.dof[i];
}


/********************************************************
* Operators for Point                                   *
* Note: we use a tolerance of 1e-8 for checking points  *
********************************************************/
bool NodeToNodeMap::Point::operator == ( const Point &rhs ) const
{
    // Two points are == if they share the same position (within tolerance)
    double dist = 0.0;
    for ( size_t i=0; i!=3; i++)
        dist += (rhs.pos[i] - pos[i]) * (rhs.pos[i] - pos[i]);
    if ( dist < 1e-16 )     // check the square of the distance (faster without sqrt)
        return true;
    return false;
}
bool NodeToNodeMap::Point::operator != ( const Point &rhs ) const
{
    return !operator==(rhs);
}
bool NodeToNodeMap::Point::operator < ( const Point &rhs ) const
{
    // Sort the points based on the x value, y value, then z-value
    for ( size_t i=0; i!=3; i++) {
        if ( (pos[i]-rhs.pos[i]) < -1e-8 ) { return true;  }
        if ( (pos[i]-rhs.pos[i]) >  1e-8 ) { return false; }
    }
    return false;
}
bool NodeToNodeMap::Point::operator <= ( const Point &rhs ) const
{
    return operator==(rhs) || operator<(rhs);
}
bool NodeToNodeMap::Point::operator >= ( const Point &rhs ) const
{
    return !operator<(rhs);
}
bool NodeToNodeMap::Point::operator > ( const Point &rhs ) const
{
    return !operator<=(rhs);
}


}
}

