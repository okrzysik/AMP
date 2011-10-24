
#include "Map3to1to3.h"
#include "Map3to1to3Parameters.h"



namespace AMP {
namespace Operator {

bool Map3to1to3::continueAsynchronousConstruction ( const boost::shared_ptr < OperatorParameters > &params )
{
    boost::shared_ptr <Map3to1to3Parameters>  p = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( params );
    AMP_ASSERT ( p );

    int size = d_MapComm.getSize();

    if ( p->d_ConstructionPhase == 0 ) {
        p->d_MasterValue = p->d_IsMaster ? 1 : 0;

        reserveRequests ( size );

        std::vector<MPI_Request>::iterator  curReq = beginRequests ();
        for (int i=0; i<size; i++) { 
            *curReq = d_MapComm.Isend( (int*) &(p->d_MasterValue), 1, i, d_SendTag );
            curReq++;
        }
        p->d_ConstructionPhase++;
        return false;

    } else if ( p->d_ConstructionPhase == 1 ) {

        p->d_NumToSend = 0;
        d_SendToProc.resize ( size );
        for (int i=0; i<size; i++) { 
            int  procDom;
            int length = d_MapComm.probe(i,d_RecvTag)/sizeof(int);
            AMP_ASSERT(length==1);
            d_MapComm.recv( (int*) &procDom, length, i, false, d_RecvTag );
            if (( procDom == 0 ) && p->d_IsMaster ) {
                d_SendToProc[i] = 1;
                p->d_NumToSend++;
            } else if (( procDom == 1 ) && !p->d_IsMaster ) {
                d_SendToProc[i] = 1;
                p->d_NumToSend++;
            } else {
                d_SendToProc[i] = 0;
            }
        }
      p->d_ConstructionPhase++;
      return false;
    } else {
      waitForAllRequests ();
      reserveRequests ( p->d_NumToSend );
    }
    return true;
}


Map3to1to3::Map3to1to3 ( const boost::shared_ptr<OperatorParameters> & params )
    : AsyncMapOperator ( params )
{
    boost::shared_ptr <Map3to1to3Parameters>  p = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( params );
    AMP_ASSERT ( p );


    d_FirstApply = true;
    d_BeginNode = p->d_SetBegin;
    d_EndNode = p->d_SetEnd;
    d_MapComm = p->d_MapComm;

    p->d_ConstructionPhase = 0;

    int DofsPerObj = p->d_db->getInteger ( "DOFsPerObject" );
    AMP_ASSERT ( ( DofsPerObj == 1 ) || ( DofsPerObj == 3 ) );

    if ( DofsPerObj == 1 ) {
        d_MapVariable = AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::NodalScalarVariable ( p->d_db->getString ( "VariableName" ) , d_MeshAdapter ) );
    }
    if ( DofsPerObj == 3 ) {
        d_MapVariable = AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::Nodal3VectorVariable ( p->d_db->getString ( "VariableName" ) , d_MeshAdapter ) );
    }

    if ( p->d_IsMaster ) {
        d_SendTag = p->d_ToSlaveCommTag;
        d_RecvTag = p->d_ToMasterCommTag;
    } else {
        d_RecvTag = p->d_ToSlaveCommTag;
        d_SendTag = p->d_ToMasterCommTag;
    }

    if ( p->d_AsynchronousConstructionParam > 0 )
        return;

    int size = d_MapComm.getSize();
    d_SendToProc.resize ( size );
    int  announceSide = p->d_IsMaster ? 1 : 0;
    std::vector<int>  whichSide ( size );
    d_MapComm.allGather( announceSide, getBufferToAvoidDebugVectorCrashing( whichSide ) );

    size_t numToSend = 0;
    for ( size_t i = 0 ; i != whichSide.size() ; i++ ) {
        if ( p->d_IsMaster && whichSide[i] == 0 ) {
            d_SendToProc[i] = 1;
            numToSend++;
        } else if ( !p->d_IsMaster && whichSide[i] == 1 ) {
            d_SendToProc[i] = 1;
            numToSend++;
        } else {
            d_SendToProc[i] = 0;
        }
    }
    reserveRequests ( numToSend );
}


void Map3to1to3::smear ( double tolerance )
{
    std::multimap<double,double>::iterator  curPtr = d_1DMultiMap.begin();
    while ( curPtr != d_1DMultiMap.end() )
    {
      std::multimap<double,double>::iterator  curLookAhead = curPtr;
      while ( curLookAhead != d_1DMultiMap.end() )
      {
        double cp = curPtr->first;
        double cl = curLookAhead->first;
        if ( fabs ( cp - cl ) > tolerance )
          break;
        curLookAhead++;
      }
      if ( curPtr != curLookAhead )
      {
        double total , count;
        total = count = 0.;
        std::multimap<double,double>::iterator  temp = curPtr;
        while ( temp != curLookAhead )
        {
          total += temp->second;
          count += 1.0;
          temp++;
        }
        curPtr->second = total / count;
        temp = curPtr;
        temp++;
        size_t  s = d_1DMultiMap.size();
        d_1DMultiMap.erase ( temp , curLookAhead );
        s = d_1DMultiMap.size();
        s++;
      }
      curPtr++;
    }
}


void Map3to1to3::addTo1DMap ( double z , double val )
{
    d_1DMultiMap.insert ( std::pair<double,double> ( z , val ) );
}


void  Map3to1to3::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr & , const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                     AMP::LinearAlgebra::Vector::shared_ptr & , const double , const double )
{
    d_1DMultiMap.clear();
    if ( d_FirstApply )
    {
      d_FirstApply = false;
    }
    else
    {
      waitForAllRequests ();
    }

    if ( !d_1DMultiMap.data().isComputed )
    {
      buildMap ( u->subsetVectorForVariable ( getInputVariable() ) );
      d_1DMultiMap.data().isComputed = true;
    }

    d_SendBuf.resize ( 2 * d_1DMultiMap.size() );
    std::vector<double>::iterator curBuf = d_SendBuf.begin();
    std::multimap<double,double>::iterator curData = d_1DMultiMap.begin();
    while ( curData != d_1DMultiMap.end() ) 
    {
      *curBuf = curData->first;
      curBuf++;
      *curBuf = curData->second;
      curBuf++;
      curData++;
    }
    std::vector<MPI_Request>::iterator  curReq = beginRequests ();
    for ( size_t i = 0 ; i != d_SendToProc.size() ; i++ )
    {
      if ( d_SendToProc[i] > 0 )
      {
        *curReq = d_MapComm.Isend( (double*) getBufferToAvoidDebugVectorCrashing( d_SendBuf ), d_SendBuf.size(), i, d_SendTag );
        curReq++;
      }
    }
}


void  Map3to1to3::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr & , const AMP::LinearAlgebra::Vector::shared_ptr & ,
                      AMP::LinearAlgebra::Vector::shared_ptr & , const double , const double )
{
    d_1DMultiMap.data().isComputed = false;
    d_1DMultiMap.clear();

    for ( size_t i = 0 ; i != d_SendToProc.size() ; i++ )
    {
      if ( d_SendToProc[i] > 0 )
      {
        int inSize = d_MapComm.probe(i,d_RecvTag)/sizeof(double);
        std::vector<double> recvBuf ( inSize );
        d_MapComm.recv( getBufferToAvoidDebugVectorCrashing(recvBuf), inSize, i, false, d_RecvTag );
        for ( int j = 0 ; j != inSize/2 ; j++ )
        {
          addTo1DMap ( recvBuf[2*j] , recvBuf[2*j+1] );
        }
      }
    }
    smear ( 1.e-8 );
    buildReturn ( d_ResultVector );

    d_ResultVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    // double a = d_ResultVector->L1Norm();

    d_SendBuf.resize ( 0 );
}


Map3to1to3::~Map3to1to3 ()
{
    waitForAllRequests ();
}


void  Map3to1to3::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &result )
{
    d_ResultVector = result->subsetVectorForVariable ( d_MapVariable );
}


void Map3to1to3::buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr p )
{
    AMP_ERROR("buildMap should never be called for the BaseClass");
}


void Map3to1to3::buildReturn ( AMP::LinearAlgebra::Vector::shared_ptr p )
{
    AMP_ERROR("buildReturn should never be called for the BaseClass");
}


}
}

