#include "Map3to1to3.h"
#include "Map3to1to3Parameters.h"
#include "ProfilerApp.h"
#include "vectors/Variable.h"


namespace AMP {
namespace Operator {


template <class T>
static T *getPtr( std::vector<T> &x )
{
    if ( x.size() == 0 ) return NULL;
    return &x[0];
}


/********************************************************
* Constructor                                           *
********************************************************/
Map3to1to3::Map3to1to3( const AMP::shared_ptr<OperatorParameters> &params_in )
    : AsyncMapOperator( params_in )
{
    // Get the input parameters
    AMP::shared_ptr<Map3to1to3Parameters> params =
        AMP::dynamic_pointer_cast<Map3to1to3Parameters>( params_in );
    AMP_ASSERT( params );
    d_commTag = params->d_commTag;

    // Create default iterators (must be overwritten by derived class)
    d_srcIterator1 = AMP::Mesh::MeshIterator();
    d_srcIterator2 = AMP::Mesh::MeshIterator();
    d_dstIterator1 = AMP::Mesh::MeshIterator();
    d_dstIterator2 = AMP::Mesh::MeshIterator();

    // Determine which processors we will be sending/receiving data from
    // 0: No communication, 1: send/recv data
    d_own_mesh1 = std::vector<bool>( d_MapComm.getSize(), false );
    d_own_mesh2 = std::vector<bool>( d_MapComm.getSize(), false );
    std::vector<char> tmp1( d_MapComm.getSize(), false );
    std::vector<char> tmp2( d_MapComm.getSize(), false );
    d_MapComm.allGather<char>( d_mesh1.get() != NULL, &tmp1[0] );
    d_MapComm.allGather<char>( d_mesh2.get() != NULL, &tmp2[0] );
    for ( int i = 0; i < d_MapComm.getSize(); i++ ) {
        d_own_mesh1[i] = tmp1[i] != 0;
        d_own_mesh2[i] = tmp2[i] != 0;
    }
    size_t numToSend = 0;
    if ( d_mesh1.get() != NULL ) {
        for ( size_t i = 0; i < d_own_mesh2.size(); i++ ) {
            if ( d_own_mesh2[i] == 1 ) numToSend++;
        }
    }
    if ( d_mesh2.get() != NULL ) {
        for ( size_t i = 0; i < d_own_mesh1.size(); i++ ) {
            if ( d_own_mesh1[i] == 1 ) numToSend++;
        }
    }
    reserveRequests( numToSend );
}


/********************************************************
* De-constructor                                        *
********************************************************/
Map3to1to3::~Map3to1to3() { waitForAllRequests(); }


/********************************************************
* Function to add points to the map                     *
********************************************************/
void Map3to1to3::addTo1DMap( std::multimap<double, double> &map, double z, double val )
{
    map.insert( std::pair<double, double>( z, val ) );
}
void Map3to1to3::addTo1DMap( std::multimap<double, double> &map,
                             const std::vector<double> &z,
                             const std::vector<double> &val )
{
    for ( size_t i = 0; i < z.size(); i++ ) map.insert( std::pair<double, double>( z[i], val[i] ) );
}


/********************************************************
* Function to start the communication                   *
********************************************************/
void Map3to1to3::applyStart( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                             AMP::LinearAlgebra::Vector::shared_ptr )
{
    const double tol = 1e-8;
    PROFILE_START( "applyStart" );

    // Subset the vector (we only need to deal with the locally owned portion)
    PROFILE_START( "subset" );
    AMP::LinearAlgebra::Variable::shared_ptr var = getInputVariable();
    AMP::LinearAlgebra::VS_Comm commSelector( AMP_MPI( AMP_COMM_SELF ) );
    AMP::LinearAlgebra::Vector::const_shared_ptr commVec =
        u->constSelect( commSelector, u->getVariable()->getName() );
    AMP::LinearAlgebra::Vector::const_shared_ptr vec = commVec->constSubsetVectorForVariable( var );
    PROFILE_STOP( "subset" );

    // Build the local maps
    PROFILE_START( "prepare data" );
    std::multimap<double, double> map1 = buildMap( vec, d_mesh1, d_srcIterator1 );
    std::multimap<double, double> map2 = buildMap( vec, d_mesh2, d_srcIterator2 );

    // Get the local 1D map coordinates
    std::multimap<double, double>::const_iterator iterator;
    double z_last = -1e100;
    std::vector<double> z1;
    for ( iterator = map1.begin(); iterator != map1.end(); ++iterator ) {
        if ( fabs( iterator->first - z_last ) > tol ) {
            z_last = iterator->first;
            z1.push_back( z_last );
        }
    }
    z_last = -1e100;
    std::vector<double> z2;
    for ( iterator = map2.begin(); iterator != map2.end(); ++iterator ) {
        if ( fabs( iterator->first - z_last ) > tol ) {
            z_last = iterator->first;
            z2.push_back( z_last );
        }
    }

    // Create the send buffers and sum the local data
    d_SendBuf1.resize( 0 ); // Reset the entries
    d_SendBuf1.resize( z1.size() );
    for ( iterator = map1.begin(); iterator != map1.end(); ++iterator ) {
        double z  = iterator->first;
        size_t i1 = std::min( AMP::Utilities::findfirst( z1, z ), z1.size() - 1 );
        size_t i2 = std::max( i1, (size_t) 1 ) - 1;
        size_t i3 = std::min( i1 + 1, z1.size() - 1 );
        size_t i  = 0;
        if ( fabs( z - z1[i1] ) < tol )
            i = i1;
        else if ( fabs( z - z1[i2] ) < tol )
            i = i2;
        else if ( fabs( z - z1[i3] ) < tol )
            i = i3;
        else
            AMP_ERROR( "Internal error" );
        d_SendBuf1[i].N++;
        d_SendBuf1[i].z = z1[i];
        d_SendBuf1[i].sum += iterator->second;
    }
    d_SendBuf2.resize( 0 ); // Reset the entries
    d_SendBuf2.resize( z2.size() );
    for ( iterator = map2.begin(); iterator != map2.end(); ++iterator ) {
        double z  = iterator->first;
        size_t i1 = std::min( AMP::Utilities::findfirst( z2, z ), z2.size() - 1 );
        size_t i2 = std::max( i1, (size_t) 1 ) - 1;
        size_t i3 = std::min( i1 + 1, z2.size() - 1 );
        size_t i  = 0;
        if ( fabs( z - z2[i1] ) < tol )
            i = i1;
        else if ( fabs( z - z2[i2] ) < tol )
            i = i2;
        else if ( fabs( z - z2[i3] ) < tol )
            i = i3;
        else
            AMP_ERROR( "Internal error" );
        d_SendBuf2[i].N++;
        d_SendBuf2[i].z = z2[i];
        d_SendBuf2[i].sum += iterator->second;
    }
    PROFILE_STOP( "prepare data" );


    // Send the data
    size_t myRank                             = (size_t) d_MapComm.getRank();
    std::vector<MPI_Request>::iterator curReq = beginRequests();
    if ( d_mesh1.get() != NULL ) {
        for ( size_t i = 0; i < d_own_mesh2.size(); i++ ) {
            if ( i == myRank ) continue; // Don't communicate local data
            if ( d_own_mesh2[i] ) {
                *curReq =
                    d_MapComm.Isend( getPtr( d_SendBuf1 ), d_SendBuf1.size(), i, d_commTag + 0 );
                ++curReq;
            }
        }
    }
    if ( d_mesh2.get() != NULL ) {
        for ( size_t i = 0; i < d_own_mesh1.size(); i++ ) {
            if ( i == myRank ) continue; // Don't communicate local data
            if ( d_own_mesh1[i] ) {
                *curReq =
                    d_MapComm.Isend( getPtr( d_SendBuf2 ), d_SendBuf2.size(), i, d_commTag + 1 );
                ++curReq;
            }
        }
    }
    PROFILE_STOP( "applyStart" );
}


/********************************************************
* Function to finish the communication and perform the  *
* interpolation                                         *
********************************************************/
void Map3to1to3::applyFinish( AMP::LinearAlgebra::Vector::const_shared_ptr,
                              AMP::LinearAlgebra::Vector::shared_ptr )
{
    PROFILE_START( "applyFinish" );

    // Recieve the data and create the maps
    size_t myRank = (size_t) d_MapComm.getRank();
    std::map<double, std::pair<int, double>> map1;
    std::map<double, std::pair<int, double>> map2;
    std::vector<comm_data> recvBuf;
    if ( d_mesh1.get() != NULL ) {
        // First get any local data
        if ( d_own_mesh2[myRank] ) unpackBuffer( d_SendBuf2, map1 );
        // Recieve all remote data
        for ( size_t i = 0; i < d_own_mesh2.size(); i++ ) {
            if ( i == myRank ) continue; // We already copied the local data
            if ( d_own_mesh2[i] ) {
                // Get the recieved data
                int inSize = d_MapComm.probe( i, d_commTag + 1 ) / sizeof( comm_data );
                recvBuf.resize( inSize );
                d_MapComm.recv( getPtr( recvBuf ), inSize, i, false, d_commTag + 1 );
                // Add it to the map
                unpackBuffer( recvBuf, map1 );
            }
        }
    }
    if ( d_mesh2.get() != NULL ) {
        // First get any local data
        if ( d_own_mesh1[myRank] ) unpackBuffer( d_SendBuf1, map2 );
        // Recieve all remote data
        for ( size_t i = 0; i < d_own_mesh1.size(); i++ ) {
            if ( i == myRank ) continue; // We already copied the local data
            if ( d_own_mesh1[i] ) {
                // Get the recieved data
                int inSize = d_MapComm.probe( i, d_commTag + 0 ) / sizeof( comm_data );
                recvBuf.resize( inSize );
                d_MapComm.recv( getPtr( recvBuf ), inSize, i, false, d_commTag + 0 );
                // Add it to the map
                unpackBuffer( recvBuf, map2 );
            }
        }
    }

    // Smear the data to create the final map
    std::map<double, std::pair<int, double>>::iterator iterator;
    std::map<double, double> final_map1;
    for ( iterator = map1.begin(); iterator != map1.end(); ++iterator ) {
        double sum = iterator->second.second;
        double N   = iterator->second.first;
        std::pair<double, double> tmp( iterator->first, sum / N );
        final_map1.insert( tmp );
    }
    std::map<double, double> final_map2;
    for ( iterator = map2.begin(); iterator != map2.end(); ++iterator ) {
        double sum = iterator->second.second;
        double N   = iterator->second.first;
        std::pair<double, double> tmp( iterator->first, sum / N );
        final_map2.insert( tmp );
    }

    // Build the return vector
    if ( d_mesh1.get() != NULL ) buildReturn( d_ResultVector, d_mesh1, d_dstIterator1, final_map1 );
    if ( d_mesh2.get() != NULL ) buildReturn( d_ResultVector, d_mesh2, d_dstIterator2, final_map2 );

    // Apply make consistent
    PROFILE_START( "makeConsistent" );
    d_ResultVector->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    PROFILE_STOP( "makeConsistent" );

    PROFILE_STOP( "applyFinish" );
}


/********************************************************
* Function to unpack a recv buffer                      *
********************************************************/
void Map3to1to3::unpackBuffer( const std::vector<comm_data> &buffer,
                               std::map<double, std::pair<int, double>> &map )
{
    const double tol = 1e-8;
    std::map<double, std::pair<int, double>>::iterator iterator, it1, it2, it3;
    for ( size_t j = 0; j < buffer.size(); j++ ) {
        iterator = map.end();
        if ( !map.empty() ) {
            it1 = map.lower_bound( buffer[j].z );
            if ( it1 == map.end() ) {
                --it1;
            }
            it2 = it1;
            it3 = it1;
            if ( it1 != map.begin() ) --it1;
            ++it3;
            if ( it3 == map.end() ) it3 = it2;
            if ( fabs( it1->first - buffer[j].z ) < tol )
                iterator = it1;
            else if ( fabs( it2->first - buffer[j].z ) < tol )
                iterator = it2;
            else if ( fabs( it3->first - buffer[j].z ) < tol )
                iterator = it3;
        }
        if ( iterator == map.end() ) {
            std::pair<int, double> tmp( buffer[j].N, buffer[j].sum );
            map.insert( std::pair<double, std::pair<int, double>>( buffer[j].z, tmp ) );
        }
        else {
            iterator->second.first += buffer[j].N;
            iterator->second.second += buffer[j].sum;
        }
    }
}


/********************************************************
* Other functions                                       *
********************************************************/
void Map3to1to3::setVector( AMP::LinearAlgebra::Vector::shared_ptr result )
{
    d_ResultVector = subsetOutputVector( result );
    AMP_ASSERT( d_ResultVector );
}


std::multimap<double, double> Map3to1to3::buildMap( AMP::LinearAlgebra::Vector::const_shared_ptr,
                                                    const AMP::Mesh::Mesh::shared_ptr,
                                                    const AMP::Mesh::MeshIterator & )
{
    AMP_ERROR( "buildMap should never be called for the BaseClass" );
    return std::multimap<double, double>();
}


void Map3to1to3::buildReturn( AMP::LinearAlgebra::Vector::shared_ptr,
                              const AMP::Mesh::Mesh::shared_ptr,
                              const AMP::Mesh::MeshIterator &,
                              const std::map<double, double> & )
{
    AMP_ERROR( "buildReturn should never be called for the BaseClass" );
}
}
}
