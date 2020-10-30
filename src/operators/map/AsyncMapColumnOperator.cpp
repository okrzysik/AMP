#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/utils/Database.h"

#include "ProfilerApp.h"

namespace AMP {
namespace Operator {


size_t globalMapTagOffset = 0; // Initialize the global map tag offset


AsyncMapColumnOperator::AsyncMapColumnOperator( const std::shared_ptr<OperatorParameters> &params )
    : AsynchronousColumnOperator( params )
{
}


void AsyncMapColumnOperator::setVector( AMP::LinearAlgebra::Vector::shared_ptr p )
{
    d_OutputVector = p;
    for ( auto &elem : d_Operators )
        std::dynamic_pointer_cast<AsyncMapOperator>( elem )->setVector( d_OutputVector );
}


void AsyncMapColumnOperator::append( std::shared_ptr<Operator> op )
{
    std::shared_ptr<AsyncMapColumnOperator> mapColumn =
        std::dynamic_pointer_cast<AsyncMapColumnOperator>( op );
    if ( mapColumn ) {
        auto curOp = mapColumn.get()->d_Operators.begin();
        while ( curOp != mapColumn.get()->d_Operators.end() ) {
            append( *curOp );
            ++curOp;
        }
    } else {
        std::shared_ptr<AsyncMapOperator> mapOp = std::dynamic_pointer_cast<AsyncMapOperator>( op );
        AMP_INSIST( mapOp, "Attempt to add a non-AsyncMapOperator to a AsyncMapColumnOperator" );
        AsynchronousColumnOperator::append( mapOp );
    }
}


void AsyncMapColumnOperator::apply( AMP::LinearAlgebra::Vector::const_shared_ptr u,
                                    AMP::LinearAlgebra::Vector::shared_ptr f )
{
    PROFILE_START( "apply" );
    this->applyStart( u, f );
    this->applyFinish( u, f );
    if ( requiresMakeConsistentSet() ) {
        AMP_ASSERT( d_OutputVector.get() != nullptr );
        d_OutputVector->makeConsistent(
            AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    }
    PROFILE_STOP( "apply" );
}


bool AsyncMapColumnOperator::requiresMakeConsistentSet()
{
    bool test = false;
    for ( auto &elem : d_Operators )
        test = test ||
               std::dynamic_pointer_cast<AsyncMapOperator>( elem )->requiresMakeConsistentSet();
    return test;
}


/********************************************************************
 * Function to copy a key from database 1 to database 2              *
 * If the key is an array of size N it will only copy the ith value. *
 ********************************************************************/
template<class TYPE>
static inline void putEntry( std::shared_ptr<AMP::Database> &database1,
                             std::shared_ptr<AMP::Database> &database2,
                             const std::string &key,
                             size_t N,
                             size_t i )
{
    auto data = database1->getVector<TYPE>( key );
    if ( N == data.size() )
        database2->putScalar( key, data[i] );
    else
        database2->putVector( key, data );
}
static void copyKey( std::shared_ptr<AMP::Database> &database1,
                     std::shared_ptr<AMP::Database> &database2,
                     std::string key,
                     int N,
                     int i )
{
    if ( database1->isDatabase( key ) ) {
        // Copy the database
        AMP_ERROR( "Not programmed for databases yet" );
    } else if ( database1->isType<bool>( key ) ) {
        // Copy a bool
        putEntry<bool>( database1, database2, key, N, i );
    } else if ( database1->isType<int>( key ) ) {
        // Copy a int
        putEntry<int>( database1, database2, key, N, i );
    } else if ( database1->isType<float>( key ) ) {
        // Copy a float
        putEntry<float>( database1, database2, key, N, i );
    } else if ( database1->isType<double>( key ) ) {
        // Copy a double
        putEntry<double>( database1, database2, key, N, i );
    } else if ( database1->isType<std::complex<double>>( key ) ) {
        // Copy a std::complex<double>
        putEntry<std::complex<double>>( database1, database2, key, N, i );
    } else if ( database1->isType<std::string>( key ) ) {
        // Copy a std::string
        putEntry<std::string>( database1, database2, key, N, i );
    } else {
        AMP_ERROR( "Unknown key type" );
    }
}


/************************************************************
 * Function to create the databases for the individual maps  *
 ************************************************************/
std::vector<std::shared_ptr<AMP::Database>>
AsyncMapColumnOperator::createDatabases( std::shared_ptr<AMP::Database> database1 )
{
    AMP_INSIST( database1->keyExists( "N_maps" ), "N_maps must exist in input database" );
    int N_maps = database1->getScalar<int>( "N_maps" );
    // Create the basic databases for each mesh
    std::vector<std::shared_ptr<AMP::Database>> meshDatabases;
    meshDatabases.reserve( N_maps );
    for ( int i = 0; i < N_maps; i++ ) {
        // Create a new database from the existing database
        std::shared_ptr<AMP::Database> database2( new AMP::Database( "MapDatabase" ) );
        std::vector<std::string> keys = database1->getAllKeys();
        for ( auto &key : keys ) {
            if ( key.compare( "N_maps" ) == 0 ) {
                // These keys are used by the builder and should not be in the sub database
            } else {
                // We need to copy the key
                copyKey( database1, database2, key, N_maps, i );
            }
        }
        meshDatabases.push_back( database2 );
    }
    return meshDatabases;
}
} // namespace Operator
} // namespace AMP
