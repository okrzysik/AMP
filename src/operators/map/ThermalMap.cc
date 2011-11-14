

#include <set>

#include "operators/map/ThermalMap.h"
#include "operators/map/ThermalMapParameters.h"


double AMP::Operator::ThermalMap::Point::_precision;

namespace AMP {
namespace Operator {
/*
  bool ThermalMap::validMapType ( const std::string &t )
  {
    if ( t == "Thermal" )
      return true;
    return false;
  }


  ThermalMap::~ThermalMap ()
  {
  }

  void ThermalMap::sendSurface ( params_shared_ptr p )
  {
    clearRequests ();
    int size = d_MapComm.getSize();

    // Get a list of the IDs on the surface
    std::set<size_t>  idsOnMySurface;
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  cur = p->d_SetBegin;
    while ( cur != p->d_SetEnd )
    {
      idsOnMySurface.insert ( cur->globalID() );
      cur++;
    }

    // Compute lists of ids and displacements to send
    p->d_ids.resize ( idsOnMySurface.size() );
    p->d_disps.resize ( idsOnMySurface.size() * 3 );

    cur = p->d_SetBegin;
    size_t curId = 0;
    while ( cur != p->d_SetEnd )
    {
      p->d_ids[curId] = cur->globalID();
      p->d_disps[3*curId + 0] = cur->x();
      p->d_disps[3*curId + 1] = cur->y();
      p->d_disps[3*curId + 2] = cur->z();
      cur++;
      curId++;
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

  void ThermalMap::recvSurface ( params_shared_ptr p )
  {
    int size = d_MapComm.getSize();
    std::multiset<Point>  surfacePts;
    for ( int i = 0 ; i != size ; i++ )
    {
      std::vector <int>     inIds;
      std::vector <double>  inDisps;
      int vecSize = d_MapComm.probe(i,d_RecvTag)/sizeof(int);
      int vecSize3 = d_MapComm.probe(i,2*d_RecvTag)/sizeof(double);
      AMP_ASSERT(vecSize3==3*vecSize);
      inIds.resize ( vecSize );
      inDisps.resize ( vecSize3 );
      d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inIds), vecSize, i, false, d_RecvTag );
      d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(inDisps), vecSize3, i, false, 2*d_RecvTag );
      for ( int j = 0 ; j != vecSize ; j++ )
      {
        Point temp;
        temp._pos[0] = inDisps [ 3*j + 0 ];
        temp._pos[1] = inDisps [ 3*j + 1 ];
        temp._pos[2] = inDisps [ 3*j + 2 ];

        bool found = false;
        std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
        while ( iter != surfacePts.upper_bound ( temp ) )
        {
          if ( (int) iter->_id == inIds[j] )
          {
            iter->_procs.push_back ( i );
            found = true;
          }
          iter++;
        }
        if ( !found )
        {
          temp._id = inIds[j];
          temp._procs.push_back ( i );
          surfacePts.insert ( temp );
        }
      }
    }

    // Calculate who receives data from me..
    std::map<size_t , CommInfo>             myIdToDestination;
    std::map<size_t , std::list<size_t> >   procToLocalIdList;
    std::map<size_t , size_t>              &remoteToLocalId = p->d_RemoteToLocalId;

    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  cur = p->d_SetBegin;
    size_t  totToSend = 0;
    while ( cur != p->d_SetEnd )
    {
      Point t;
      t._pos[0] = cur->x();
      t._pos[1] = cur->y();
      t._pos[2] = cur->z();
      std::multiset<Point>::iterator whichPt = surfacePts.lower_bound ( t );
      AMP_ASSERT ( whichPt != surfacePts.end() );
      AMP_ASSERT ( whichPt != surfacePts.upper_bound ( t ) );
      AMP_ASSERT ( *whichPt == t );
      remoteToLocalId [ whichPt->_id ] = cur->globalID();
      CommInfo &data = myIdToDestination [ cur->globalID() ];
      data._remId = whichPt->_id;
      data._procs.insert ( data._procs.begin() , whichPt->_procs.begin() , whichPt->_procs.end() );
      std::list<int>::iterator curProc = data._procs.begin();
      while ( curProc != data._procs.end() )
      {
        procToLocalIdList[*curProc].push_back ( cur->globalID() );
        curProc++;
        totToSend++;
      }
      cur++;
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
      std::list<size_t>::iterator  localToSend = procToLocalIdList[i].begin();
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
  }

  void ThermalMap::sendOrder ( params_shared_ptr p )
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

  void ThermalMap::recvOrder ( params_shared_ptr p )
  {
    int size = d_MapComm.getSize();
    for ( int i = 0 ; i != size ; i++ )
    {
      int *buf = getBufferToAvoidDebugVectorCrashing ( d_MySurfaceIndicesRecv );
      d_MapComm.recv( (int*) buf+d_MySurfaceDisplsR[i], d_MySurfaceCountsR[i], i, false, d_RecvTag );
    }
  }

  void ThermalMap::buildSendRecvList ( params_shared_ptr p )
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

  void ThermalMap::finalizeCommunication ( params_shared_ptr p )
  {
    waitForAllRequests ();
    clearRequests ();
    reserveRequests ( 2*p->d_NumPartners );
  }

  bool ThermalMap::continueAsynchronousConstruction ( const boost::shared_ptr < AMP::Operator::OperatorParameters > &params )
  {
    params_shared_ptr p = boost::dynamic_pointer_cast<ThermalMapParameters> ( params );
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

  ThermalMap::Point::Point ()
  {
  }

  ThermalMap::Point::Point ( const Point &rhs )
  {
    for ( size_t i = 0 ; i != 3 ; i++ )
      _pos[i] = rhs._pos[i];
    _id = rhs._id;
    _procs.insert ( _procs.begin() , rhs._procs.begin() , rhs._procs.end() );
  }

  bool ThermalMap::Point::operator == ( const Point &rhs ) const
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

  bool ThermalMap::Point::operator < ( const Point &rhs ) const
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


  ThermalMap::ThermalMap ( const boost::shared_ptr<AMP::Operator::OperatorParameters> & params )
       : AMP::Operator::AsyncMapOperator ( params )
  {
    // Cast the params appropriately
    d_OutputVector = AMP::LinearAlgebra::Vector::shared_ptr ();
    AMP_ASSERT ( params );
    ThermalMapParameters &Params =
                        *(boost::dynamic_pointer_cast<ThermalMapParameters> ( params ) );


    // Set class members
    d_MapComm = Params.d_MapComm;
    int DofsPerObj = Params.d_db->getInteger ( "DOFsPerObject" );
    AMP_ASSERT ( ( DofsPerObj == 1 ) || ( DofsPerObj == 3 ) );

    if ( DofsPerObj == 1 )
    {
      d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::NodalScalarVariable ( Params.d_db->getString ( "VariableName" ) , d_MeshAdapter ) );
    }
    if ( DofsPerObj == 3 )
    {
      d_inpVariable = AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::Nodal3VectorVariable ( Params.d_db->getString ( "VariableName" ) , d_MeshAdapter ) );
    }

    int commSize = d_MapComm.getSize();
    d_MySurfaceCountsS.resize ( commSize );
    d_MySurfaceDisplsS.resize ( commSize );
    d_MySurfaceCountsR.resize ( commSize );
    d_MySurfaceDisplsR.resize ( commSize );
    Point::_precision = 1.e-8;

    // If this needs to be built asynchronously, exit
    if ( Params.d_AsynchronousConstructionParam > 0 )
    {
      return;
    }

    d_SendTag = d_RecvTag = Params.d_ToMasterCommTag;

    // Get a list of the IDs on the surface
    std::set<size_t>  idsOnMySurface;
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  cur = Params.d_SetBegin;
    while ( cur != Params.d_SetEnd )
    {
      idsOnMySurface.insert ( cur->globalID() );
      cur++;
    }

    // Compute lists of ids and displacements to send
    std::vector<int>     ids ( idsOnMySurface.size() );
    std::vector<double>  disps ( idsOnMySurface.size() * 3 );

    cur = Params.d_SetBegin;
    size_t curId = 0;
    while ( cur != Params.d_SetEnd )
    {
      ids[curId] = cur->globalID();
      disps[3*curId + 0] = cur->x();
      disps[3*curId + 1] = cur->y();
      disps[3*curId + 2] = cur->z();
      cur++;
      curId++;
    }

    // Send out my surface to everyone else
    int  numMyVals = curId;
    d_MapComm.allGather( numMyVals, getBufferToAvoidDebugVectorCrashing(d_MySurfaceCountsS) );
    d_MySurfaceDisplsS[0] = 0;
    for (int i=1; i!=commSize; i++)
    {
      d_MySurfaceDisplsS[i] = d_MySurfaceDisplsS[i-1] + d_MySurfaceCountsS[i-1];
    }
    std::vector<int>    otherIds ( d_MySurfaceDisplsS[commSize-1] +
                                   d_MySurfaceCountsS[commSize-1] );
    std::vector<double> otherDisps ( otherIds.size() * 3 );

    d_MapComm.allGather( getBufferToAvoidDebugVectorCrashing( ids ), 
                         ids.size(), 
                         getBufferToAvoidDebugVectorCrashing( otherIds ),
                         getBufferToAvoidDebugVectorCrashing( d_MySurfaceCountsS ),
                         getBufferToAvoidDebugVectorCrashing( d_MySurfaceDisplsS ),
                         true );

    for (int i=0; i!=commSize; i++)
    {
      d_MySurfaceCountsS[i] *= 3;
      d_MySurfaceDisplsS[i] *= 3;
    }
    d_MapComm.allGather( getBufferToAvoidDebugVectorCrashing ( disps ),
                         disps.size(),
                         getBufferToAvoidDebugVectorCrashing ( otherDisps ),
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceCountsS ),
                         getBufferToAvoidDebugVectorCrashing ( d_MySurfaceDisplsS ),
                         true );

    // Parse the surfaces into a sorted list
    std::multiset<Point>  surfacePts;
    for (int i=0; i!=commSize; i++)
    {
      int numToRead = d_MySurfaceCountsS[i]/3;
      int offsetIDs = d_MySurfaceDisplsS[i]/3;
      int offsetDisps = d_MySurfaceDisplsS[i];
      for (int j=0; j!=numToRead; j++ )
      {
        Point temp;
        temp._pos[0] = otherDisps [ offsetDisps + 3*j + 0 ];
        temp._pos[1] = otherDisps [ offsetDisps + 3*j + 1 ];
        temp._pos[2] = otherDisps [ offsetDisps + 3*j + 2 ];

        bool found = false;
        std::multiset<Point>::iterator  iter = surfacePts.lower_bound ( temp );
        while ( iter != surfacePts.upper_bound ( temp ) )
        {
          if ( (int) iter->_id == otherIds[offsetIDs + j] )
          {
            iter->_procs.push_back ( i );
            found = true;
          }
          iter++;
        }
        if ( !found )
        {
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

    cur = Params.d_SetBegin;
    size_t  totToSend = 0;
    while ( cur != Params.d_SetEnd )
    {
      Point t;
      t._pos[0] = cur->x();
      t._pos[1] = cur->y();
      t._pos[2] = cur->z();
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
      while ( curProc != data._procs.end() )
      {
        procToLocalIdList[*curProc].push_back ( cur->globalID() );
        curProc++;
        totToSend++;
      }
      cur++;
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

  void ThermalMap::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
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

  void  ThermalMap::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &p )
  {
    d_OutputVector = p->subsetVectorForVariable ( d_inpVariable );
    AMP_INSIST ( d_OutputVector , "setVector received bogus stuff" );
  }

  void ThermalMap::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr &  ,
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
*/
}
}

