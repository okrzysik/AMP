#include "operators/map/NodeToNodeMap.h"
#include "operators/map/NodeToNodeMapParameters.h"
#include "ampmesh/MeshElement.h"
#include "discretization/NodalVariable.h"

#include <set>


namespace AMP {
namespace Operator {


double AMP::Operator::NodeToNodeMap::Point::_precision;


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
    int dim = d_Mesh->getDim();
    AMP_INSIST(dim==3,"Node to Node map has only been tested in 3d (see Point)");
    d_MapComm = Params.d_MapComm;
    int DofsPerObj = Params.d_db->getInteger ( "DOFsPerObject" );
    AMP_ASSERT ( ( DofsPerObj == 1 ) || ( DofsPerObj == 3 ) );

    // Create a nodal variable 
    AMP::LinearAlgebra::Variable::shared_ptr variable( new AMP::Discretization::NodalVariable(DofsPerObj,"VariableName") );    

    int commSize = d_MapComm.getSize();
    d_MySurfaceCountsS.resize ( commSize );
    d_MySurfaceDisplsS.resize ( commSize );
    d_MySurfaceCountsR.resize ( commSize );
    d_MySurfaceDisplsR.resize ( commSize );
    Point::_precision = 1.e-8;

    // If this needs to be built asynchronously, exit
    if ( Params.d_AsynchronousConstructionParam > 0 )
        return;

    d_SendTag = d_RecvTag = Params.d_ToMasterCommTag;

    // Compute lists of ids and displacements to send
    AMP::Mesh::MeshIterator  cur = Params.d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator  end = cur.end();
    std::vector<AMP::Mesh::MeshElementID>  ids( cur.size() );
    std::vector<double>  disps( cur.size()*dim );
    for (size_t i=0; i<ids.size(); i++) {
        ids[i] = cur->globalID();
        std::vector<double> pos = cur->coord();
        AMP_ASSERT((int)pos.size()==dim);
        for (int j=0; j<dim; j++)
            disps[dim*i+j] = pos[j];
        ++cur;
    }

    // Send out my surface to everyone else
    int  numMyVals = ids.size();
    d_MapComm.allGather( numMyVals, getBufferToAvoidDebugVectorCrashing(d_MySurfaceCountsS) );
    d_MySurfaceDisplsS[0] = 0;
    for (int i=1; i!=commSize; i++)
        d_MySurfaceDisplsS[i] = d_MySurfaceDisplsS[i-1] + d_MySurfaceCountsS[i-1];
    int recieveCount = d_MySurfaceDisplsS[commSize-1] + d_MySurfaceCountsS[commSize-1];
    std::vector<AMP::Mesh::MeshElementID>  otherIds( recieveCount );
    std::vector<double> otherDisps( recieveCount*dim );

    d_MapComm.allGather( getBufferToAvoidDebugVectorCrashing( ids ), 
                         ids.size(), 
                         getBufferToAvoidDebugVectorCrashing( otherIds ),
                         getBufferToAvoidDebugVectorCrashing( d_MySurfaceCountsS ),
                         getBufferToAvoidDebugVectorCrashing( d_MySurfaceDisplsS ),
                         true );

    for (int i=0; i!=commSize; i++) {
        d_MySurfaceCountsS[i] *= dim;
        d_MySurfaceDisplsS[i] *= dim;
    }
    d_MapComm.allGather( getBufferToAvoidDebugVectorCrashing ( disps ),
                         disps.size(),
                         getBufferToAvoidDebugVectorCrashing ( otherDisps ),
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceCountsS ),
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceDisplsS ),
                         true );

    // Parse the surfaces into a sorted list
    std::multiset<Point>  surfacePts;
    for (int i=0; i!=commSize; i++) {
        int numToRead = d_MySurfaceCountsS[i]/dim;
        int offsetIDs = d_MySurfaceDisplsS[i]/dim;
        int offsetDisps = d_MySurfaceDisplsS[i];
        for (int j=0; j!=numToRead; j++ ) {
            Point temp;
            for (int k=0; k<dim; k++)
                temp._pos[k] = otherDisps[ offsetDisps + dim*j + k ];
            bool found = false;
            std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
            while ( iter != surfacePts.upper_bound ( temp ) ) {
                if ( iter->_id == otherIds[offsetIDs + j] ) {
                    iter->_procs.push_back ( i );
                    found = true;
                }
                iter++;
            }
            if ( !found ) {
                temp._id = otherIds[offsetIDs + j];
                temp._procs.push_back ( i );
                surfacePts.insert ( temp );
            }
        }
    }

    // Calculate who receives data from me..
    std::map<size_t , CommInfo>             myIdToDestination;
    std::map<size_t , std::list<size_t> >   procToLocalIdList;
    std::map<size_t , size_t>               remoteToLocalId;

    cur = cur.begin();
    size_t  totToSend = 0;
    while ( cur != end ) {
        Point t;
        std::vector<double> pos = cur->coord();
        for (int j=0; j<dim; j++)
            t._pos[j] = pos[j];
        std::multiset<Point>::iterator whichPt = surfacePts.lower_bound ( t );
        AMP_INSIST ( whichPt != surfacePts.end() , "Node-to-node map not aligned correctly" );
        AMP_INSIST ( whichPt != surfacePts.upper_bound ( t ) , "Node-to-node map not aligned correctly" );
        if ( whichPt->_id == cur->globalID() ) whichPt++;
        AMP_INSIST ( whichPt != surfacePts.end() , "Node-to-node map not aligned correctly" );
        AMP_INSIST ( *whichPt == t , "Node-to-node map not aligned correctly" );

        remoteToLocalId [ whichPt->_id ] = cur->globalID();
        CommInfo &data = myIdToDestination [ cur->globalID() ];
        data._remId = whichPt->_id;
        data._procs.insert ( data._procs.begin() , whichPt->_procs.begin() , whichPt->_procs.end() );
        std::list<int>::iterator curProc = data._procs.begin();
        while ( curProc != data._procs.end() ) {
            procToLocalIdList[*curProc].push_back ( cur->globalID() );
            curProc++;
            totToSend++;
        }
        ++cur;
    }

    // Compute the send/recv vectors for the all to all communication
    d_MySurfaceIndicesSend.resize ( DofsPerObj * totToSend );
    d_SendBuffer.resize ( DofsPerObj * totToSend );
    d_RecvBuffer.resize ( DofsPerObj * totToSend );
    int curOffset = 0;
    for (int i = 0; i!=commSize; i++)
    {
      d_MySurfaceDisplsR[i] = d_MySurfaceDisplsS[i] = curOffset;
      std::list<size_t>::iterator  localToSend = procToLocalIdList[i].begin();
      while ( localToSend != procToLocalIdList[i].end() )
      {
        for (int j=0; j!=DofsPerObj; j++)
        {
          d_MySurfaceIndicesSend[curOffset++] = *localToSend;
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

}


/********************************************************
* De-constructor                                        *
********************************************************/
NodeToNodeMap::~NodeToNodeMap ()
{
}



bool NodeToNodeMap::validMapType ( const std::string &t )
{
    if ( t == "NodeToNode" )
        return true;
    return false;
}


void NodeToNodeMap::sendSurface ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
    clearRequests ();
    int size = d_MapComm.getSize();
    int dim = d_Mesh->getDim();

    // Compute lists of ids and displacements to send
    AMP::Mesh::MeshIterator  cur = p->d_BoundaryNodeIterator.begin();
    AMP::Mesh::MeshIterator  end = p->d_BoundaryNodeIterator.end();
    p->d_ids.resize( cur.size() );
    p->d_disps.resize( cur.size()*dim );

    cur =  p->d_BoundaryNodeIterator.begin();
    for (size_t i=0; i<p->d_ids.size(); i++) {
        p->d_ids[i] = cur->globalID();
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
        *curReq = d_MapComm.Isend( getBufferToAvoidDebugVectorCrashing(p->d_ids), p->d_ids.size(), i, d_SendTag );
        curReq++;
        *curReq = d_MapComm.Isend( getBufferToAvoidDebugVectorCrashing(p->d_disps), p->d_disps.size(), i, 2*d_SendTag );
        curReq++;
    }
}


void NodeToNodeMap::recvSurface ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
    int size = d_MapComm.getSize();
    int dim = d_Mesh->getDim();
    AMP_INSIST(dim==3,"NodeToNodeMap needs to be fixed for dimensions other than 3");
    std::multiset<Point>  surfacePts;
    for (int i=0; i!=size; i++) {
        std::vector <AMP::Mesh::MeshElementID>     inIds;
        std::vector <double>  inDisps;
        int vecSize = d_MapComm.probe(i,d_RecvTag)/sizeof(int);
        int vecSize3 = d_MapComm.probe(i,2*d_RecvTag)/sizeof(double);
        AMP_ASSERT(vecSize3==3*vecSize);
        inIds.resize ( vecSize );
        inDisps.resize ( vecSize3 );
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inIds), vecSize, i, false, d_RecvTag );
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inDisps), vecSize3, i, false, 2*d_RecvTag );
        for ( int j = 0 ; j != vecSize ; j++ ) {
            Point temp;
            temp._pos[0] = inDisps [ 3*j + 0 ];
            temp._pos[1] = inDisps [ 3*j + 1 ];
            temp._pos[2] = inDisps [ 3*j + 2 ];

            bool found = false;
            std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
            while ( iter != surfacePts.upper_bound ( temp ) ) {
                if ( iter->_id == inIds[j] ) {
                    iter->_procs.push_back ( i );
                    found = true;
                }
                iter++;
            }
            if ( !found ) {
                temp._id = inIds[j];
                temp._procs.push_back ( i );
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
            t._pos[j] = pos[j];
        std::multiset<Point>::iterator whichPt = surfacePts.lower_bound ( t );
        AMP_ASSERT ( whichPt != surfacePts.end() );
        AMP_ASSERT ( whichPt != surfacePts.upper_bound ( t ) );
        AMP_ASSERT ( *whichPt == t );
        //remoteToLocalId [ whichPt->_id ] = cur->globalID();
        CommInfo &data = myIdToDestination [ cur->globalID() ];
        data._remId = whichPt->_id;
        data._procs.insert ( data._procs.begin() , whichPt->_procs.begin() , whichPt->_procs.end() );
        std::list<int>::iterator curProc = data._procs.begin();
        while ( curProc != data._procs.end() ) {
            procToLocalIdList[*curProc].push_back ( cur->globalID() );
            curProc++;
            totToSend++;
        }
        ++cur;
    }
AMP_ERROR("Not finished converting");
/*
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
    waitForAllRequests ();
    clearRequests ();
    reserveRequests ( d_MySurfaceCountsS.size() );
    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    for ( size_t i = 0 ; i != d_MySurfaceCountsS.size() ; i++ )
    {
      int *buf = getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesSend );
      *curReq = d_MapComm.Isend( (int*) buf + d_MySurfaceDisplsS[i],d_MySurfaceCountsS[i], i, d_SendTag );
      curReq++;
    }
}


void NodeToNodeMap::recvOrder ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
    int size = d_MapComm.getSize();
    for ( int i = 0 ; i != size ; i++ )
    {
      int *buf = getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesRecv );
      d_MapComm.recv( (int*) buf+d_MySurfaceDisplsR[i], d_MySurfaceCountsR[i], i, false, d_RecvTag );
    }
}


void NodeToNodeMap::buildSendRecvList ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
    std::map<size_t , size_t>              &remoteToLocalId = p->d_RemoteToLocalId;
    int DofsPerObj = p->d_db->getInteger ( "DOFsPerObject" );

    AMP::Mesh::DOFMap::shared_ptr  varDof = d_MeshAdapter->getDOFMap ( d_inpVariable );
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
}


void NodeToNodeMap::finalizeCommunication ( boost::shared_ptr<NodeToNodeMapParameters> p )
{
    waitForAllRequests ();
    clearRequests ();
    reserveRequests ( 2*p->d_NumPartners );
}


bool NodeToNodeMap::continueAsynchronousConstruction ( const boost::shared_ptr < AMP::Operator::OperatorParameters > &params )
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
}


NodeToNodeMap::Point::Point ()
{
}


NodeToNodeMap::Point::Point ( const Point &rhs )
{
    for ( size_t i = 0 ; i != 3 ; i++ )
      _pos[i] = rhs._pos[i];
    _id = rhs._id;
    _procs.insert ( _procs.begin() , rhs._procs.begin() , rhs._procs.end() );
}


bool NodeToNodeMap::Point::operator == ( const Point &rhs ) const
{
    double dist = 0.0;
    for ( size_t i = 0 ; i != 3 ; i++ )
    {
      dist += (rhs._pos[i] - _pos[i]) * (rhs._pos[i] - _pos[i]);
    }
    if ( sqrt ( dist ) < _precision )
      return true;
    return false;
}


bool NodeToNodeMap::Point::operator < ( const Point &rhs ) const
{
    if ( !operator == ( rhs ) )
    {
      for ( size_t i = 0 ; i != 3 ; i++ )
      {
        if ( ( _pos[i] - rhs._pos[i] ) > _precision ) return false;
        if ( ( rhs._pos[i] - _pos[i] ) > _precision ) return true;
      }
    }
    return false;
}



void NodeToNodeMap::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
                              const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                                    AMP::LinearAlgebra::Vector::shared_ptr & ,
                              const double ,
                              const double )
{
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
}


}
}

