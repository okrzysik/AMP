#include "AMP/ampmesh/loadBalance/loadBalanceSimulator.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MultiMesh.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>


// inline function to convert two size_t to double, divide, round the result, and return a size_t
// integer
static inline size_t divide_double( size_t x, size_t y )
{
    return (size_t) floor( ( ( (double) x ) / ( (double) y ) ) + 0.5 );
}


namespace AMP {
namespace Mesh {


/************************************************************
 * Constructors                                              *
 ************************************************************/
loadBalanceSimulator::loadBalanceSimulator()
{
    d_N_elements = 0;
    d_decomp     = 0;
    d_min        = 0;
    d_max        = 0;
    d_max_ranks  = 0;
    cache_valid  = false;
}
loadBalanceSimulator::loadBalanceSimulator( AMP::shared_ptr<MeshParameters> params,
                                            const std::vector<int> &ranks,
                                            size_t N_elements )
{
    // Get required values from the parameters
    AMP_ASSERT( !ranks.empty() );
    AMP::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT( database != nullptr );
    AMP_INSIST( database->keyExists( "MeshType" ), "MeshType must exist in input database" );
    AMP_INSIST( database->keyExists( "MeshName" ), "MeshName must exist in input database" );
    d_name = database->getString( "MeshName" );
    d_type = database->getString( "MeshType" );
    // Simulate the load process
    d_decomp    = 0;
    cache_valid = false;
    if ( d_type == std::string( "Multimesh" ) ) {
        loadBalanceSimulator rhs = MultiMesh::simulateBuildMesh( params, ranks );
        d_name                   = rhs.d_name;
        d_type                   = rhs.d_type;
        d_N_elements             = rhs.d_N_elements;
        d_max_ranks              = rhs.d_max_ranks;
        d_params                 = rhs.d_params;
        d_ranks                  = rhs.d_ranks;
        d_submeshes              = rhs.d_submeshes;
        d_decomp                 = rhs.d_decomp;
        cache_valid              = rhs.cache_valid;
        if ( d_ranks.size() == 1 ) {
            d_min       = d_N_elements;
            d_max       = d_N_elements;
            cache_valid = true;
        }
    } else {
        d_params = params;
        d_ranks  = ranks;
        d_submeshes.resize( 0 );
        if ( N_elements == 0 )
            d_N_elements = Mesh::estimateMeshSize( params );
        else
            d_N_elements = N_elements;
        d_min       = divide_double( d_N_elements, d_ranks.size() );
        d_max       = divide_double( d_N_elements, d_ranks.size() );
        d_max_ranks = Mesh::maxProcs( params );
        cache_valid = true;
    }
    if ( d_ranks.size() > d_max_ranks )
        d_ranks.resize( d_max_ranks );
}
loadBalanceSimulator::loadBalanceSimulator( AMP::shared_ptr<MeshParameters> params,
                                            const std::vector<int> &ranks,
                                            const std::vector<loadBalanceSimulator> &submeshes,
                                            int decomp )
{
    AMP::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT( database != nullptr );
    AMP_INSIST( database->keyExists( "MeshType" ), "MeshType must exist in input database" );
    AMP_INSIST( database->keyExists( "MeshName" ), "MeshName must exist in input database" );
    d_name       = database->getString( "MeshName" );
    d_type       = database->getString( "MeshType" );
    d_params     = params;
    d_ranks      = ranks;
    d_submeshes  = submeshes;
    d_decomp     = decomp;
    cache_valid  = false;
    d_N_elements = 0;
    d_max_ranks  = ~static_cast<size_t>( 0 );
    for ( auto &elem : d_submeshes )
        d_N_elements += elem.d_N_elements;
    if ( d_ranks.size() == 1 ) {
        d_min       = d_N_elements;
        d_max       = d_N_elements;
        cache_valid = true;
    }
}
loadBalanceSimulator::loadBalanceSimulator( const loadBalanceSimulator &rhs )
    : d_name( rhs.d_name ),
      d_type( rhs.d_type ),
      d_ranks( rhs.d_ranks ),
      d_submeshes( rhs.d_submeshes )
{
    d_N_elements = rhs.d_N_elements;
    d_max_ranks  = rhs.d_max_ranks;
    d_params     = rhs.d_params;
    d_decomp     = rhs.d_decomp;
    cache_valid  = rhs.cache_valid;
    d_min        = rhs.d_min;
    d_max        = rhs.d_max;
}


/************************************************************
 * Function to add a processor                               *
 ************************************************************/
bool loadBalanceSimulator::addProc( int rank )
{
    bool added  = false;
    cache_valid = false;
    if ( d_submeshes.empty() ) {
        if ( d_ranks.size() < d_max_ranks ) {
            d_ranks.push_back( rank );
            added = true;
        }
    } else if ( d_type == std::string( "Multimesh" ) ) {
        added = MultiMesh::addProcSimulation( *this, d_submeshes, rank, d_decomp );
        if ( added )
            d_ranks.push_back( rank );
    } else {
        AMP_ERROR( "Not finished" );
    }
    return added;
}


/************************************************************
 * Function to return the min, max, and avg # of elements    *
 ************************************************************/
size_t loadBalanceSimulator::min()
{
    if ( !cache_valid )
        updateCache();
    return d_min;
}
size_t loadBalanceSimulator::max()
{
    if ( !cache_valid )
        updateCache();
    return d_max;
}
size_t loadBalanceSimulator::avg() { return divide_double( d_N_elements, d_ranks.size() ); }


/************************************************************
 * Function to print the load balance                        *
 ************************************************************/
void loadBalanceSimulator::print( unsigned char detail, unsigned char indent_N )
{
    int N_procs = 0;
    for ( auto &elem : d_ranks )
        N_procs = std::max( N_procs, elem + 1 );
    char indent[257];
    memset( indent, 0, 257 );
    memset( indent, 0x20, indent_N );
    if ( detail & 0x1 ) {
        // Print the global load balance
        std::vector<size_t> N_elements( N_procs, 0 );
        countElements( *this, N_elements );
        std::cout << indent << "Rank, N_elements:" << std::endl;
        int N_line = 16;
        for ( int i = 0; i < ( N_procs + N_line - 1 ) / N_line; i++ ) {
            std::cout << indent;
            for ( int j = i * N_line; j < std::min( ( i + 1 ) * N_line, N_procs ); j++ )
                std::cout << std::setw( 8 ) << j;
            std::cout << std::endl;
            std::cout << indent;
            for ( int j = i * N_line; j < std::min( ( i + 1 ) * N_line, N_procs ); j++ )
                std::cout << std::setw( 8 ) << N_elements[j];
            std::cout << std::endl << std::endl;
        }
    }
    if ( detail & 0x2 ) {
        std::cout << indent << d_name << ": " << d_ranks.size() << std::endl;
        for ( auto &elem : d_submeshes )
            elem.print( 2, indent_N + 3 );
    }
}


/************************************************************
 * Misc. functions                                           *
 ************************************************************/
void loadBalanceSimulator::changeRanks( const std::vector<int> &ranks )
{
    if ( d_submeshes.empty() ) {
        d_ranks = ranks;
        updateCache();
    } else {
        AMP_INSIST( ranks.size() == d_ranks.size(), "Cannot change ranks to different size" );
        for ( size_t i = 0; i < ranks.size(); i++ )
            d_ranks[i] = ranks[i];
    }
}
void loadBalanceSimulator::countElements( const loadBalanceSimulator &mesh,
                                          std::vector<size_t> &N_elements )
{
    if ( mesh.d_submeshes.empty() ) {
        for ( size_t i = 0; i < mesh.d_ranks.size(); i++ )
            N_elements[mesh.d_ranks[i]] += mesh.d_N_elements / mesh.d_ranks.size();
    } else {
        for ( auto &elem : mesh.d_submeshes )
            countElements( elem, N_elements );
    }
}
void loadBalanceSimulator::updateCache()
{
    if ( d_submeshes.empty() ) {
        d_min = divide_double( d_N_elements, d_ranks.size() );
        d_max = divide_double( d_N_elements, d_ranks.size() );
        return;
    }
    if ( d_decomp == 0 ) {
        // General case
        int N_procs = 0;
        for ( auto &elem : d_ranks )
            N_procs = std::max( N_procs, elem + 1 );
        std::vector<size_t> N_elements( N_procs, 0 );
        countElements( *this, N_elements );
        d_min = d_N_elements;
        d_max = 0;
        for ( auto &elem : d_ranks ) {
            d_min = std::min( d_min, N_elements[elem] );
            d_max = std::max( d_max, N_elements[elem] );
        }
    } else if ( d_decomp == 1 ) {
        // Special case where no two submeshes share a processor
        d_min = d_N_elements;
        d_max = 0;
        for ( auto &elem : d_submeshes ) {
            if ( elem.cache_valid ) {
                if ( elem.d_min < d_min )
                    d_min = elem.d_min;
                if ( elem.d_max > d_max )
                    d_max = elem.d_max;
            } else {
                d_min = std::min( d_min, elem.min() );
                d_max = std::max( d_max, elem.max() );
            }
        }
    } else {
        AMP_ERROR( "Unkown value for d_decomp" );
    }
    cache_valid = true;
}
} // namespace Mesh
} // namespace AMP
