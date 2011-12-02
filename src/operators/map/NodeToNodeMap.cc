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
    int commSize = d_MapComm.getSize();
    DofsPerObj = Params.d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj<=3,"Node to Node map only works for <= 3 DOFs per node (see Point)");

    // Create a nodal variable and DOFManager (this should be moved out of here)
    d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr( new AMP::Discretization::NodalVariable(DofsPerObj,"VariableName") );
    if ( d_mesh1.get() != NULL )
        d_DOFManager1 = AMP::Discretization::DOFManager::shared_ptr( new AMP::Discretization::simpleDOFManager(d_mesh1,AMP::Mesh::Vertex,0,DofsPerObj) );
    if ( d_mesh2.get() != NULL )
        d_DOFManager2 = AMP::Discretization::DOFManager::shared_ptr( new AMP::Discretization::simpleDOFManager(d_mesh2,AMP::Mesh::Vertex,0,DofsPerObj) );

    // For each mesh, get the list of points owned by the current processor
    d_ownedPointsMesh1 = std::vector<Point>();
    d_ownedPointsMesh2 = std::vector<Point>();
    if ( d_mesh1.get() != NULL )
        d_ownedPointsMesh1 = createOwnedPoints( d_mesh1->getIDsetIterator(AMP::Mesh::Vertex,Params.d_BoundaryID1,0), d_DOFManager1 );
    if ( d_mesh2.get() != NULL )
        d_ownedPointsMesh2 = createOwnedPoints( d_mesh2->getIDsetIterator(AMP::Mesh::Vertex,Params.d_BoundaryID2,0), d_DOFManager2 );

    // Send the list of points on mesh1 to all processors
    int send_cnt = (int) d_ownedPointsMesh1.size();
    std::vector<int> recv_cnt(commSize,0);
    std::vector<int> recv_disp(commSize,0);
    d_MapComm.allGather(send_cnt,&recv_cnt[0]);
    for (int i=1; i<commSize; i++)
        recv_disp[i] = recv_disp[i-1]+recv_cnt[i-1];
    int N_recv_tot = recv_disp[commSize-1] + recv_cnt[commSize-1];
    std::vector<Point> surfacePts = std::vector<Point>(N_recv_tot);
    d_MapComm.allGather( getPtr(d_ownedPointsMesh1), send_cnt, &surfacePts[0], &recv_cnt[0], &recv_disp[0], true );

    // Sort the points for fast searching
    AMP::Utilities::quicksort(surfacePts);

    // Find the points in mesh1 that align with the points owned by the current processor on mesh2
    std::vector<Point> mesh1Point(d_ownedPointsMesh2.size());
    for (size_t i=0; i<d_ownedPointsMesh2.size(); i++) {
        size_t index = AMP::Utilities::findfirst(surfacePts,d_ownedPointsMesh2[i]);
        AMP_ASSERT(surfacePts[index]==d_ownedPointsMesh2[i]);
        mesh1Point[i] = surfacePts[index];
    }
    surfacePts.clear();


    AMP_ERROR("Not finished");

/*
    // Determine which processors own which meshes 
    // This will be needed to determine which processors are involved in which communications
    d_owner = std::vector<unsigned char>(commSize,0);
    unsigned char local_owner = 0;
    if ( d_mesh1.get()!=NULL )
        local_owner += 1;
    if ( d_mesh2.get()!=NULL )
        local_owner += 2;
    d_MapComm.allGather(local_owner,&d_owner[0]);


    // Parse the surfaces into a sorted list
    std::multiset<Point>  surfacePts;
    for (int i=0; i!=commSize; i++) {
        int numToRead = d_MySurfaceCountsS[i]/dim;
        int offsetIDs = d_MySurfaceDisplsS[i]/dim;
        int offsetDisps = d_MySurfaceDisplsS[i];
        for (int j=0; j!=numToRead; j++ ) {
            Point temp;
            for (int k=0; k<dim; k++)
                temp.pos[k] = otherDisps[ offsetDisps + dim*j + k ];
            bool found = false;
            std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
            while ( iter != surfacePts.upper_bound ( temp ) ) {
                if ( iter->id == otherIds[offsetIDs + j] ) {
                    iter->procs.push_back ( i );
                    found = true;
                }
                iter++;
            }
            if ( !found ) {
                temp.id = otherIds[offsetIDs + j];
                temp.procs.push_back ( i );
                surfacePts.insert ( temp );
            }
        }
    }

    // Calculate who receives data from me..
    std::map<AMP::Mesh::MeshElementID,CommInfo>            myIdToDestination;
    std::map< int, std::list<AMP::Mesh::MeshElementID> >   procToLocalIdList;
    //std::map<size_t , size_t>               remoteToLocalId;

    cur = cur.begin();
    size_t  totToSend = 0;
    while ( cur != end ) {
        Point t;
        std::vector<double> pos = cur->coord();
        for (int j=0; j<dim; j++)
            t.pos[j] = pos[j];
        std::multiset<Point>::iterator whichPt = surfacePts.lower_bound ( t );
        AMP_INSIST ( whichPt != surfacePts.end() , "Node-to-node map not aligned correctly" );
        AMP_INSIST ( whichPt != surfacePts.upper_bound ( t ) , "Node-to-node map not aligned correctly" );
        if ( whichPt->id == cur->globalID() ) whichPt++;
        AMP_INSIST ( whichPt != surfacePts.end() , "Node-to-node map not aligned correctly" );
        AMP_INSIST ( *whichPt == t , "Node-to-node map not aligned correctly" );

        //remoteToLocalId [ whichPt->id ] = cur->globalID();
        CommInfo &data = myIdToDestination [ cur->globalID() ];
        data._remId = whichPt->id;
        data.procs.insert ( data.procs.begin() , whichPt->procs.begin() , whichPt->procs.end() );
        std::list<int>::iterator curProc = data.procs.begin();
        while ( curProc != data.procs.end() ) {
            procToLocalIdList[*curProc].push_back ( cur->globalID() );
            curProc++;
            totToSend++;
        }
        ++cur;
    }

    // Compute the list of DOFs and the displacements for all nodes on the surface id
    d_MyDOFsSend.resize ( DofsPerObj * totToSend );
    d_SendBuffer.resize ( DofsPerObj * totToSend );
    d_RecvBuffer.resize ( DofsPerObj * totToSend );
    int curOffset = 0;
    for (int i = 0; i!=commSize; i++)
    {
      d_MySurfaceDisplsR[i] = d_MySurfaceDisplsS[i] = curOffset;
      std::list<AMP::Mesh::MeshElementID>::iterator  localToSend = procToLocalIdList[i].begin();
      while ( localToSend != procToLocalIdList[i].end() ) {
        for (int j=0; j!=DofsPerObj; j++)
        {
          d_MyDOFsSend[curOffset++] = *localToSend;
        }
        localToSend++;
      }
      d_MySurfaceCountsR[i] = d_MySurfaceCountsS[i] = curOffset - d_MySurfaceDisplsS[i];
    }

    d_MySurfaceIndicesRecv.resize ( DofsPerObj * totToSend );
    d_MapComm.allToAll ( getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesSend ) ,
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceCountsS ) ,
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceDisplsS ) ,
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesRecv ) ,
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceCountsR ) ,
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceDisplsR ) ,
                         true );

    AMP::Mesh::DOFMap::shared_ptr  varDof = d_MeshAdapter->getDOFMap ( d_inpVariable );
    for ( size_t i = 0 ; i != d_MySurfaceIndicesSend.size() ; i++ )
    {
      AMP_ASSERT ( remoteToLocalId.find ( d_MySurfaceIndicesRecv[i] ) != remoteToLocalId.end() );
      d_MySurfaceIndicesRecv[i] = varDof->getGlobalID ( remoteToLocalId[d_MySurfaceIndicesRecv[i]] , i % DofsPerObj );
      d_MySurfaceIndicesSend[i] = varDof->getGlobalID ( d_MySurfaceIndicesSend[i] , i % DofsPerObj );
    }

    size_t numPartners = 0;
    for ( size_t i = 0 ; i != d_MySurfaceCountsS.size() ; i++ )
    {
      if ( d_MySurfaceCountsS[i] > 0 )
      {
        numPartners++;
      }
    }
    reserveRequests ( 2*numPartners );
*/
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



void NodeToNodeMap::sendSurface ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*    clearRequests ();
    int size = d_MapComm.getSize();
    int dim = d_Mesh->getDim();

    // Compute lists of ids and displacements to send
    AMP::Mesh::MeshIterator  cur = p->d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator  end = p->d_BoundaryNodeIterator.end();
    p->dids.resize( cur.size() );
    p->d_disps.resize( cur.size()*dim );

    cur =  p->d_BoundaryNodeIterator.begin();
    for (size_t i=0; i<p->dids.size(); i++) {
        p->dids[i] = cur->globalID();
        std::vector<double> pos = cur->coord();
        for (int j=0; j<dim; j++)
            p->d_disps[dim*i+j] = pos[j];
        ++cur;
    }

    if ( p->d_IsMaster )
    {
      d_SendTag = p->d_ToSlaveCommTag;
      d_RecvTag = p->d_ToMasterCommTag;
    }
    else
    {
      d_SendTag = p->d_ToMasterCommTag;
      d_RecvTag = p->d_ToSlaveCommTag;
    }

    reserveRequests ( 2*size );
    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    for ( int i = 0 ; i != size ; i++ )
    {
        *curReq = d_MapComm.Isend( getBufferToAvoidDebugVectorCrashing(p->dids), p->dids.size(), i, d_SendTag );
        curReq++;
        *curReq = d_MapComm.Isend( getBufferToAvoidDebugVectorCrashing(p->d_disps), p->d_disps.size(), i, 2*d_SendTag );
        curReq++;
    }
*/
}


void NodeToNodeMap::recvSurface ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*
    int size = d_MapComm.getSize();
    int dim = d_Mesh->getDim();
    AMP_INSIST(dim==3,"NodeToNodeMap needs to be fixed for dimensions other than 3");
    std::multiset<Point>  surfacePts;
    for (int i=0; i!=size; i++) {
        std::vector <AMP::Mesh::MeshElementID>     inIds;
        std::vector <double>  inDisps;
        int vecSize = d_MapComm.probe(i,d_RecvTag)/sizeof(AMP::Mesh::MeshElementID);
        int vecSizeX = d_MapComm.probe(i,2*d_RecvTag)/sizeof(double);
        AMP_ASSERT(vecSizeX==dim*vecSize);
        inIds.resize ( vecSize );
        inDisps.resize ( vecSizeX );
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inIds), vecSize, i, false, d_RecvTag );
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inDisps), vecSizeX, i, false, 2*d_RecvTag );
        for (int j=0; j!=vecSize; j++) {
            Point temp;
            for (int k=0; k<dim; k++)
                temp.pos[k] = inDisps[dim*j+k];

            bool found = false;
            std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
            while ( iter != surfacePts.upper_bound ( temp ) ) {
                if ( iter->id == inIds[j] ) {
                    iter->procs.push_back ( i );
                    found = true;
                }
                iter++;
            }
            if ( !found ) {
                temp.id = inIds[j];
                temp.procs.push_back ( i );
                surfacePts.insert ( temp );
            }
        }
    }

    // Calculate who receives data from me..
    std::map< AMP::Mesh::MeshElementID, CommInfo >   myIdToDestination;
    std::map< int, std::list<AMP::Mesh::MeshElementID> >   procToLocalIdList;
    //std::map<size_t , size_t>              &remoteToLocalId = p->d_RemoteToLocalId;
    AMP::Mesh::MeshIterator  cur = p->d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator  end = p->d_BoundaryNodeIterator.end();
    size_t  totToSend = 0;
    while ( cur != end ) {
        Point t;
        std::vector<double> pos = cur->coord();
        for (int j=0; j<dim; j++)
            t.pos[j] = pos[j];
        std::multiset<Point>::iterator whichPt = surfacePts.lower_bound(t);
        AMP_ASSERT ( whichPt != surfacePts.end() );             // Check that iterator points to a valid point
        AMP_ASSERT ( whichPt != surfacePts.upper_bound(t) );    // Check that the point was found
        AMP_ASSERT ( *whichPt == t );                           // Check that the point matches the found point
        //remoteToLocalId [ whichPt->id ] = cur->globalID();
        CommInfo &data = myIdToDestination [ cur->globalID() ];
        data._remId = whichPt->id;
        data.procs.insert ( data.procs.begin() , whichPt->procs.begin() , whichPt->procs.end() );
        std::list<int>::iterator curProc = data.procs.begin();
        while ( curProc != data.procs.end() ) {
            procToLocalIdList[*curProc].push_back ( cur->globalID() );
            curProc++;
            totToSend++;
        }
        ++cur;
    }

    int DofsPerObj = p->d_db->getInteger ( "DOFsPerObject" );

    // Compute the send/recv vectors for the all to all communication
    d_MySurfaceIndicesSend.resize ( DofsPerObj * totToSend );
    d_MySurfaceIndicesRecv.resize ( DofsPerObj * totToSend );
    d_SendBuffer.resize ( DofsPerObj * totToSend );
    d_RecvBuffer.resize ( DofsPerObj * totToSend );
    int curOffset = 0;
    for (int i=0; i!=size; i++ )
    {
      d_MySurfaceDisplsR[i] = d_MySurfaceDisplsS[i] = curOffset;
      std::list<AMP::Mesh::MeshElementID>::iterator  localToSend = procToLocalIdList[i].begin();
      while ( localToSend != procToLocalIdList[i].end() )
      {
        for (int j=0; j!=DofsPerObj; j++ )
        {
          d_MySurfaceIndicesSend[curOffset++] = *localToSend;
        }
        localToSend++;
      }
      d_MySurfaceCountsR[i] = d_MySurfaceCountsS[i] = curOffset - d_MySurfaceDisplsS[i];
    }
*/
}


void NodeToNodeMap::sendOrder ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*
    waitForAllRequests ();
    clearRequests ();
    reserveRequests ( d_MySurfaceCountsS.size() );
    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    for ( size_t i = 0 ; i != d_MySurfaceCountsS.size() ; i++ )
    {
        AMP::Mesh::MeshElementID* buf = NULL;
        if ( d_MySurfaceIndicesSend.size()>0 )
            buf = &d_MySurfaceIndicesSend[d_MySurfaceDisplsS[i]];
        *curReq = d_MapComm.Isend( (int*) buf + d_MySurfaceDisplsS[i], d_MySurfaceCountsS[i], i, d_SendTag );
        curReq++;
    }
*/
}


void NodeToNodeMap::recvOrder ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*
    int size = d_MapComm.getSize();
    for ( int i = 0 ; i != size ; i++ )
    {
        AMP::Mesh::MeshElementID* buf = NULL;
        if ( d_MySurfaceIndicesRecv.size()>0 )
            buf = &d_MySurfaceIndicesRecv[d_MySurfaceDisplsR[i]];
        d_MapComm.recv( (int*) buf+d_MySurfaceDisplsR[i], d_MySurfaceCountsR[i], i, false, d_RecvTag );
    }
*/
}


void NodeToNodeMap::buildSendRecvList ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*
    std::map<size_t , size_t>              &remoteToLocalId = p->d_RemoteToLocalId;
    int DofsPerObj = p->d_db->getInteger ( "DOFsPerObject" );

    for ( size_t i = 0 ; i != d_MySurfaceIndicesSend.size() ; i++ )
    {
      AMP_ASSERT ( remoteToLocalId.find ( d_MySurfaceIndicesRecv[i] ) != remoteToLocalId.end() );
      d_MySurfaceIndicesRecv[i] = varDof->getGlobalID ( remoteToLocalId[d_MySurfaceIndicesRecv[i]] , i % DofsPerObj );
      d_MySurfaceIndicesSend[i] = varDof->getGlobalID ( d_MySurfaceIndicesSend[i] , i % DofsPerObj );
    }

    p->d_NumPartners = 0;
    for ( size_t i = 0 ; i != d_MySurfaceCountsS.size() ; i++ )
    {
      if ( d_MySurfaceCountsS[i] > 0 )
      {
        p->d_NumPartners++;
      }
    }
*/
}


void NodeToNodeMap::finalizeCommunication ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
/*
    waitForAllRequests ();
    clearRequests ();
    reserveRequests ( 2*p->d_NumPartners );
*/
}


/*bool NodeToNodeMap::continueAsynchronousConstruction ( const boost::shared_ptr < AMP::Operator::OperatorParameters > &params )
{
    boost::shared_ptr<NodeToNodeMapParameters> p = boost::dynamic_pointer_cast<NodeToNodeMapParameters> ( params );
    AMP_ASSERT ( p );

    switch ( p->d_ConstructionPhase )
    {
      case  0:  sendSurface ( p );
                break;

      case  1:  recvSurface ( p );
                break;

      case  2:  sendOrder ( p );
                break;

      case  3:  recvOrder ( p );
                break;

      case  4:  buildSendRecvList ( p );
                break;

      case  5:  finalizeCommunication ( p );
                break;

      case  6:  return true; // Done!

      default:  AMP_INSIST ( 1 , "Bad construction phase" );
    }
    p->d_ConstructionPhase++;
    return false;
}*/



void NodeToNodeMap::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
                              const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                                    AMP::LinearAlgebra::Vector::shared_ptr & ,
                              const double ,
                              const double )
{
AMP_ERROR("Not converted yet");
/*
    AMP::LinearAlgebra::Vector::shared_ptr   curPhysics = u->subsetVectorForVariable ( d_inpVariable );
    AMP_INSIST ( curPhysics , "apply received bogus stuff" );

    curPhysics->getValuesByGlobalID ( d_MySurfaceIndicesSend.size() ,
                                      getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesSend ) ,
                                      getBufferToAvoidDebugVectorCrashing ( d_SendBuffer ) );

    std::vector<MPI_Request>::iterator  curReq = beginRequests();
    double *buf1 = getBufferToAvoidDebugVectorCrashing ( d_SendBuffer );
    double *buf2 = getBufferToAvoidDebugVectorCrashing ( d_RecvBuffer );
    for ( size_t i = 0 ; i != d_MySurfaceCountsS.size() ; i++ )
    {
      if ( d_MySurfaceCountsS[i] > 0 )
      {
        *curReq = d_MapComm.Isend( static_cast<double*>( buf1 + d_MySurfaceDisplsS[i] ), d_MySurfaceCountsS[i], i, d_SendTag );
        curReq++;
        *curReq = d_MapComm.Irecv( static_cast<double*>( buf2 + d_MySurfaceDisplsR[i] ), d_MySurfaceCountsR[i], i, d_RecvTag );
        curReq++;
      }
    }
*/
}


void  NodeToNodeMap::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p )
{
    d_OutputVector = p->subsetVectorForVariable ( d_inpVariable );
    AMP_INSIST ( d_OutputVector , "setVector received bogus stuff" );
}


void NodeToNodeMap::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
                              const AMP::LinearAlgebra::Vector::shared_ptr &   ,
                                    AMP::LinearAlgebra::Vector::shared_ptr &r  ,
                              const double ,
                              const double )
{
AMP_ERROR("Not converted yet");
/*

    waitForAllRequests ();
//    bool  copyToOutput = false;

    AMP::LinearAlgebra::Vector::shared_ptr   curPhysics;
//    if ( r )
//    {
//      curPhysics = r->subsetVectorForVariable ( d_inpVariable );
//    }
//    else if ( d_OutputVector )
//    {
//      copyToOutput = false;
//      curPhysics = d_OutputVector;
//    }

    double *buf = getBufferToAvoidDebugVectorCrashing ( d_RecvBuffer );
//    curPhysics->setValuesByGlobalID ( d_MySurfaceIndicesRecv.size() ,
//                                      getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesRecv ) ,
//                                      buf );

//    if ( d_OutputVector && copyToOutput )
//    {
      d_OutputVector->setValuesByGlobalID ( d_MySurfaceIndicesRecv.size() ,
                                      getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesRecv ) ,
                                      buf );
//    }

    d_OutputVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
*/
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
        Point temp;
        temp.id = cur->globalID();
        std::vector<double> pos = cur->coord();
        DOFManager->getDOFs(temp.id,dofs);
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

