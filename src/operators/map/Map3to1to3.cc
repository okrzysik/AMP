#include "Map3to1to3.h"
#include "Map3to1to3Parameters.h"
#include "vectors/Variable.h"
#include "utils/ProfilerApp.h"


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
Map3to1to3::Map3to1to3 ( const boost::shared_ptr<OperatorParameters> & params_in )
    : AsyncMapOperator ( params_in )
{
    // Get the input parameters
    boost::shared_ptr <Map3to1to3Parameters>  params = boost::dynamic_pointer_cast<Map3to1to3Parameters> ( params_in );
    AMP_ASSERT ( params );
    d_commTag = params->d_commTag;
    int DofsPerObj = params->d_db->getInteger ( "DOFsPerObject" );
    AMP_INSIST(DofsPerObj==1,"Map3to1to3 is currently only designed for 1 DOF per node");

    // Create default iterators (must be overwritten by derived class)
    d_srcIterator1 = AMP::Mesh::MeshIterator();
    d_srcIterator2 = AMP::Mesh::MeshIterator();
    d_dstIterator1 = AMP::Mesh::MeshIterator();
    d_dstIterator2 = AMP::Mesh::MeshIterator();

    // Determine which processors we will be sending/receiving data from
    // 0: No communication, 1: send/recv data
    d_own_mesh1 = std::vector<bool>(d_MapComm.getSize(),false);
    d_own_mesh2 = std::vector<bool>(d_MapComm.getSize(),false);
    std::vector<char> tmp1(d_MapComm.getSize(),false);
    std::vector<char> tmp2(d_MapComm.getSize(),false);
    d_MapComm.allGather<char>( d_mesh1.get()!=NULL, &tmp1[0] );
    d_MapComm.allGather<char>( d_mesh2.get()!=NULL, &tmp2[0] );
    for (int i=0; i<d_MapComm.getSize(); i++) {
        d_own_mesh1[i] = tmp1[i]!=0;
        d_own_mesh2[i] = tmp2[i]!=0;
    }
    size_t numToSend = 0;
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh2.size(); i++) {
            if ( d_own_mesh2[i]==1 )
                numToSend++;
        }
    }
    if ( d_mesh2.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh1.size(); i++) {
            if ( d_own_mesh1[i]==1 )
                numToSend++;
        }
    }
    reserveRequests ( numToSend );
}


/********************************************************
* De-constructor                                        *
********************************************************/
Map3to1to3::~Map3to1to3 ()
{
    waitForAllRequests ();
}


/********************************************************
* Function to smear nearby points                       *
********************************************************/
void Map3to1to3::smear( std::multimap<double,double> &map, double tolerance )
{
    std::multimap<double,double>::iterator  curPtr = map.begin();
    while ( curPtr != map.end() ) {
        std::multimap<double,double>::iterator  curLookAhead = curPtr;
        while ( curLookAhead != map.end() ) {
            double cp = curPtr->first;
            double cl = curLookAhead->first;
            if ( fabs ( cp - cl ) > tolerance )
              break;
            curLookAhead++;
        }
        if ( curPtr != curLookAhead ) {
            double total , count;
            total = count = 0.;
            std::multimap<double,double>::iterator  temp = curPtr;
            while ( temp != curLookAhead ) {
                total += temp->second;
                count += 1.0;
                temp++;
            }
            curPtr->second = total / count;
            temp = curPtr;
            temp++;
            size_t  s = map.size();
            map.erase( temp , curLookAhead );
            s = map.size();
            s++;
        }
        curPtr++;
    }
}


/********************************************************
* Function to add points to the map                     *
********************************************************/
void Map3to1to3::addTo1DMap ( std::multimap<double,double> &map, double z , double val )
{
    map.insert ( std::pair<double,double> ( z , val ) );
}


/********************************************************
* Function to start the communication                   *
********************************************************/
void  Map3to1to3::applyStart ( const AMP::LinearAlgebra::Vector::shared_ptr & , const AMP::LinearAlgebra::Vector::shared_ptr &u ,
                     AMP::LinearAlgebra::Vector::shared_ptr & , const double , const double )
{
    PROFILE_START("applyStart");

    // Build the local maps
    AMP::LinearAlgebra::Vector::shared_ptr vec = subsetInputVector( u );
    std::multimap<double,double> map1 = buildMap( vec, d_mesh1, d_srcIterator1 );
    std::multimap<double,double> map2 = buildMap( vec, d_mesh2, d_srcIterator2 );

    // Create the send buffers
    d_SendBuf1.resize(2*map1.size());
    std::multimap<double,double>::iterator curData = map1.begin();
    for (size_t i=0; i<map1.size(); i++) {
        d_SendBuf1[2*i+0] = curData->first;
        d_SendBuf1[2*i+1] = curData->second;
        curData++;
    }
    d_SendBuf2.resize(2*map2.size());
    curData = map2.begin();
    for (size_t i=0; i<map2.size(); i++) {
        d_SendBuf2[2*i+0] = curData->first;
        d_SendBuf2[2*i+1] = curData->second;
        curData++;
    }

    // Send the data
    std::vector<MPI_Request>::iterator  curReq = beginRequests();
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh2.size(); i++ ) {
            if ( d_own_mesh2[i] ) {
                *curReq = d_MapComm.Isend( getPtr(d_SendBuf1), d_SendBuf1.size(), i, d_commTag+0 );
                curReq++;
            }
        }
    }
    if ( d_mesh2.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh1.size(); i++ ) {
            if ( d_own_mesh1[i] ) {
                *curReq = d_MapComm.Isend( getPtr(d_SendBuf2), d_SendBuf2.size(), i, d_commTag+1 );
                curReq++;
            }
        }
    }
    PROFILE_STOP("applyStart");
}


/********************************************************
* Function to finish the communication and perform the  *
* interpolation                                         *
********************************************************/
void  Map3to1to3::applyFinish ( const AMP::LinearAlgebra::Vector::shared_ptr & , const AMP::LinearAlgebra::Vector::shared_ptr & ,
                      AMP::LinearAlgebra::Vector::shared_ptr & , const double , const double )
{
    PROFILE_START("applyFinish");
    // Recieve the data and create the maps
    std::multimap<double,double> map1;
    std::multimap<double,double> map2;
    std::vector<double> recvBuf;
    if ( d_mesh1.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh2.size(); i++ ) {
            if ( d_own_mesh2[i] ) {
                int inSize = d_MapComm.probe(i,d_commTag+1)/sizeof(double);
                recvBuf.resize( inSize );
                d_MapComm.recv( getPtr(recvBuf), inSize, i, false, d_commTag+1 );
                for (int j=0; j<inSize/2; j++)
                    addTo1DMap( map1, recvBuf[2*j], recvBuf[2*j+1] );
            }
        }
    }
    if ( d_mesh2.get() != NULL ) {
        for (size_t i=0; i<d_own_mesh1.size(); i++ ) {
            if ( d_own_mesh1[i] ) {
                int inSize = d_MapComm.probe(i,d_commTag+0)/sizeof(double);
                recvBuf.resize( inSize );
                d_MapComm.recv( getPtr(recvBuf), inSize, i, false, d_commTag+0 );
                for (int j=0; j<inSize/2; j++)
                    addTo1DMap( map2, recvBuf[2*j], recvBuf[2*j+1] );
            }
        }
    }

    // Smear the data
    if ( d_mesh1.get() != NULL )
        smear( map1, 1.e-8 );
    if ( d_mesh2.get() != NULL )
        smear( map2, 1.e-8 );

    // Build the return vector
    if ( d_mesh1.get() != NULL )
        buildReturn( d_ResultVector, d_mesh1, d_dstIterator1, map1 );
    if ( d_mesh2.get() != NULL )
        buildReturn( d_ResultVector, d_mesh2, d_dstIterator2, map2 );

    // Apply make consistent
    d_ResultVector->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    // double a = d_ResultVector->L1Norm();
    PROFILE_STOP("applyFinish");
}



/********************************************************
* Other functions                                       *
********************************************************/
void  Map3to1to3::setVector ( AMP::LinearAlgebra::Vector::shared_ptr &result )
{
    d_ResultVector = subsetInputVector( result );
    AMP_ASSERT ( d_ResultVector );
}


std::multimap<double,double>  Map3to1to3::buildMap ( const AMP::LinearAlgebra::Vector::shared_ptr, 
    const AMP::Mesh::Mesh::shared_ptr, const AMP::Mesh::MeshIterator& )
{
    AMP_ERROR("buildMap should never be called for the BaseClass");
    return std::multimap<double,double>();
}


void Map3to1to3::buildReturn ( AMP::LinearAlgebra::Vector::shared_ptr, const AMP::Mesh::Mesh::shared_ptr, 
    const AMP::Mesh::MeshIterator&, const std::multimap<double,double>& )
{
    AMP_ERROR("buildReturn should never be called for the BaseClass");
}


}
}

