#include "AMP/operators/map/AsyncMapColumnOperator.h"
#include "AMP/operators/map/AsyncMapOperator.h"
#include "AMP/utils/Database.h"

#include "ProfilerApp.h"

namespace AMP::Operator {


size_t globalMapTagOffset = 0; // Initialize the global map tag offset


/********************************************************
 * Constructors                                          *
 ********************************************************/
AsyncMapColumnOperator::AsyncMapColumnOperator() : AsynchronousColumnOperator() {}
AsyncMapColumnOperator::AsyncMapColumnOperator( std::shared_ptr<const OperatorParameters> params )
    : AsynchronousColumnOperator( params )
{
}


void AsyncMapColumnOperator::setVector( AMP::LinearAlgebra::Vector::shared_ptr p )
{
    d_OutputVector = p;
    for ( auto &elem : d_operators )
        std::dynamic_pointer_cast<AsyncMapOperator>( elem )->setVector( d_OutputVector );
}


void AsyncMapColumnOperator::append( std::shared_ptr<Operator> op )
{
    auto mapColumn = std::dynamic_pointer_cast<AsyncMapColumnOperator>( op );
    if ( mapColumn ) {
        auto curOp = mapColumn.get()->d_operators.begin();
        while ( curOp != mapColumn.get()->d_operators.end() ) {
            append( *curOp );
            ++curOp;
        }
    } else {
        auto mapOp = std::dynamic_pointer_cast<AsyncMapOperator>( op );
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
        AMP_ASSERT( d_OutputVector );
        d_OutputVector->makeConsistent(
            AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    }
    PROFILE_STOP( "apply" );
}


bool AsyncMapColumnOperator::requiresMakeConsistentSet()
{
    bool test = false;
    for ( auto &elem : d_operators )
        test = test ||
               std::dynamic_pointer_cast<AsyncMapOperator>( elem )->requiresMakeConsistentSet();
    return test;
}


/********************************************************************
 * Function to copy a key from database 1 to database 2              *
 * If the key is an array of size N it will only copy the ith value. *
 ********************************************************************/
template<class TYPE>
static inline void putEntry( std::shared_ptr<AMP::Database> database1,
                             std::shared_ptr<AMP::Database> database2,
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
static void copyKey( std::shared_ptr<AMP::Database> database1,
                     std::shared_ptr<AMP::Database> database2,
                     std::string key,
                     int N,
                     int i )
{
    if ( database1->isDatabase( key ) ) {
        // Copy the database
        auto db = database1->getDatabase( key );
        database2->putDatabase( key, db->cloneDatabase() );
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
    int N_maps = database1->getScalar<int>( "N_maps" );
    // Create the basic databases for each mesh
    std::vector<std::shared_ptr<AMP::Database>> meshDatabases;
    meshDatabases.reserve( N_maps );
    for ( int i = 0; i < N_maps; i++ ) {
        // Create a new database from the existing database
        auto database2 = std::make_shared<AMP::Database>( "MapDatabase" );
        auto keys      = database1->getAllKeys();
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
} // namespace AMP::Operator
