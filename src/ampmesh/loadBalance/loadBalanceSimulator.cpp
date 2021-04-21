#include "AMP/ampmesh/loadBalance/loadBalanceSimulator.h"
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>


namespace AMP::Mesh {


template<class TYPE>
static inline TYPE min( const std::vector<TYPE> &x )
{
    TYPE y = std::numeric_limits<TYPE>::max();
    for ( auto i : x )
        y = std::min( y, i );
    return y;
}
template<class TYPE>
static inline TYPE max( const std::vector<TYPE> &x )
{
    TYPE y = 0;
    for ( auto i : x )
        y = std::max( y, i );
    return y;
}
template<class TYPE>
static inline TYPE sum( const std::vector<TYPE> &x )
{
    TYPE y = 0;
    for ( TYPE i : x )
        y += i;
    return y;
}


template<class TYPE>
static inline int findMin( const std::vector<TYPE> &x )
{
    int i  = -1;
    TYPE y = std::numeric_limits<TYPE>::max();
    for ( int j = 0; j < (int) x.size(); j++ ) {
        if ( x[j] < y ) {
            i = j;
            y = x[j];
        }
    }
    return i;
}
template<class TYPE>
static inline int findMax( const std::vector<TYPE> &x )
{
    int i  = -1;
    TYPE y = 0;
    for ( int j = 0; j < (int) x.size(); j++ ) {
        if ( x[j] > y ) {
            i = j;
            y = x[j];
        }
    }
    return i;
}


static inline bool allEqual( const std::vector<loadBalanceSimulator> &data )
{
    bool equal = true;
    for ( size_t i = 1; i < data.size(); i++ )
        equal = equal && data[0].getCost() == data[i].getCost();
    return equal;
}


/********************************************************
 * Divide M items with weights into N groups ( M > N )   *
 ********************************************************/
static std::vector<int> divideGroups( int N, const std::vector<double> &x )
{
    // Divide the entries into N groups minimizing the maximum cost
    AMP_ASSERT( (int) x.size() >= N );
    std::vector<std::pair<double, int>> ids( x.size() );
    for ( size_t i = 0; i < x.size(); i++ )
        ids[i] = std::make_pair( x[i], i );
    AMP::Utilities::quicksort( ids );
    std::vector<double> cost( N, 0 );
    std::vector<int> groups( x.size(), -1 );
    while ( !ids.empty() ) {
        // Find the entry with the lowest cost
        int i = findMin( cost );
        // Add the largest remaining entry to the lowest cost group
        cost[i] += ids.back().first;
        groups[ids.back().second] = i;
        ids.resize( ids.size() - 1 );
    }
    return groups;
}


/************************************************************
 * Constructors                                              *
 ************************************************************/
loadBalanceSimulator::loadBalanceSimulator()
    : d_cost( 0 ), d_maxCostRank( 0 ), d_max_procs( 0 ), d_method( 0 )
{
}
loadBalanceSimulator::loadBalanceSimulator( std::shared_ptr<AMP::Database> db )
{
    // Get required values from the parameters
    AMP_ASSERT( db );
    AMP_INSIST( db->keyExists( "MeshType" ), "MeshType must exist in input database" );
    AMP_INSIST( db->keyExists( "MeshName" ), "MeshName must exist in input database" );
    d_name    = db->getString( "MeshName" );
    auto type = db->getString( "MeshType" );
    // Simulate the load process
    if ( type == std::string( "Multimesh" ) ) {
        // Create a database for each mesh within the multimesh
        auto meshDatabases = MultiMesh::createDatabases( db );
        d_submeshes.resize( meshDatabases.size() );
        for ( size_t i = 0; i < meshDatabases.size(); i++ )
            d_submeshes[i] = loadBalanceSimulator( meshDatabases[i] );
        d_cost = 0;
        for ( size_t i = 0; i < meshDatabases.size(); i++ )
            d_cost += d_submeshes[i].d_cost;
        d_method = db->getWithDefault( "LoadBalanceMethod", 2 );
        if ( d_method == 0 ) {
            d_max_procs = 100000000;
            for ( size_t i = 0; i < meshDatabases.size(); i++ )
                d_max_procs = std::min( d_max_procs, d_submeshes[i].d_max_procs );
        } else {
            d_max_procs = 0;
            for ( size_t i = 0; i < meshDatabases.size(); i++ )
                d_max_procs += d_submeshes[i].d_max_procs;
        }
    } else {
        d_submeshes.resize( 0 );
        auto params = std::make_shared<MeshParameters>( db );
        d_max_procs = Mesh::maxProcs( params );
        d_cost      = Mesh::estimateMeshSize( params );
        d_method    = 0;
    }
    d_cost *= db->getWithDefault( "Weight", 1.0 );
    d_maxCostRank = d_cost;
    AMP_ASSERT( d_max_procs > 0 );
    AMP_ASSERT( d_cost > 0 );
}
loadBalanceSimulator::loadBalanceSimulator( double cost, int maxProc, const std::string &name )
    : d_name( name ), d_cost( cost ), d_maxCostRank( cost ), d_max_procs( maxProc ), d_method( 0 )
{
    AMP_ASSERT( d_cost > 0 );
    if ( d_max_procs == 0 )
        d_max_procs = 1000000;
}
loadBalanceSimulator::loadBalanceSimulator( const std::vector<loadBalanceSimulator> &meshes,
                                            int method,
                                            const std::string &name )
    : d_name( name ),
      d_cost( 0 ),
      d_maxCostRank( 0 ),
      d_max_procs( 0 ),
      d_method( method ),
      d_submeshes( meshes )
{
    for ( size_t i = 0; i < d_submeshes.size(); i++ )
        d_cost += d_submeshes[i].d_cost;
    d_maxCostRank = d_cost;
    if ( d_method == 0 ) {
        d_max_procs = 100000000;
        for ( size_t i = 0; i < d_submeshes.size(); i++ )
            d_max_procs = std::min( d_max_procs, d_submeshes[i].d_max_procs );
    } else {
        d_max_procs = 0;
        for ( size_t i = 0; i < d_submeshes.size(); i++ )
            d_max_procs += d_submeshes[i].d_max_procs;
    }
}


/************************************************************
 * Function to return the min, max, and avg # of elements    *
 ************************************************************/
std::vector<double> loadBalanceSimulator::getRankCost() const
{
    PROFILE_START( "getRankCost", 1 );
    int N_proc = 0;
    for ( auto rank : d_ranks )
        N_proc = std::max( N_proc, rank + 1 );
    std::vector<double> cost( N_proc, 0 );
    addRankCost( cost );
    PROFILE_STOP( "getRankCost", 1 );
    return cost;
}
void loadBalanceSimulator::addRankCost( std::vector<double> &cost ) const
{
    if ( !d_submeshes.empty() ) {
        for ( const auto &mesh : d_submeshes )
            mesh.addRankCost( cost );
    } else {
        double N = d_cost / d_ranks.size();
        for ( auto rank : d_ranks ) {
            AMP_ASSERT( rank < (int) cost.size() );
            cost[rank] += N;
        }
    }
}


/************************************************************
 * Function to print the load balance                        *
 ************************************************************/
void loadBalanceSimulator::print( uint8_t detail, uint8_t indent_N )
{
    int N_procs = 0;
    for ( auto &elem : d_ranks )
        N_procs = std::max( N_procs, elem + 1 );
    char indent[512] = { 0 };
    memset( indent, 0x20, indent_N );
    if ( detail & 0x1 ) {
        // Print the global load balance
        auto N_elements = getRankCost();
        std::cout << indent << "Rank, N_elements:" << std::endl;
        int N_line = 16;
        for ( int i = 0; i < ( N_procs + N_line - 1 ) / N_line; i++ ) {
            std::cout << indent;
            for ( int j = i * N_line; j < std::min( ( i + 1 ) * N_line, N_procs ); j++ )
                std::cout << std::setw( 8 ) << j;
            std::cout << std::endl;
            std::cout << indent;
            for ( int j = i * N_line; j < std::min( ( i + 1 ) * N_line, N_procs ); j++ )
                std::cout << std::setw( 8 ) << static_cast<int>( N_elements[j] );
            std::cout << std::endl << std::endl;
        }
    }
    if ( detail & 0x2 ) {
        printf( "%s%s (%0.0f): %i", indent, d_name.data(), d_cost, (int) d_ranks.size() );
        if ( !d_ranks.empty() ) {
            printf( " (%i", d_ranks[0] );
            for ( size_t i = 1; i < d_ranks.size() && i < 10; i++ )
                printf( ",%i", d_ranks[i] );
            if ( d_ranks.size() > 10 )
                printf( ",..." );
            printf( ")" );
        }
        printf( "\n" );
        for ( auto &elem : d_submeshes )
            elem.print( 2, indent_N + 3 );
    }
}


/************************************************************
 * Misc. functions                                           *
 ************************************************************/
void loadBalanceSimulator::setProcs( int N_procs )
{
    std::vector<int> ranks( N_procs );
    for ( int i = 0; i < N_procs; i++ )
        ranks[i] = i;
    setRanks( ranks );
}
static void setCost0( const std::vector<int> &N,
                      const std::vector<loadBalanceSimulator> &mesh,
                      std::vector<double> &cost )
{
    int cost0 = 0;
    for ( size_t i = 0; i < N.size(); i++ ) {
        if ( N[i] == 0 )
            cost0 += mesh[i].getCost();
    }
    for ( size_t i = 0; i < N.size(); i++ ) {
        if ( N[i] == 0 )
            cost[i] = cost0;
    }
}
void loadBalanceSimulator::loadBalance( int N_proc, std::vector<int> &N )
{
    PROFILE_SCOPED( timer, "setRanks-loadBalance", 2 );
    // Initialize the cost
    std::vector<double> cost( d_submeshes.size(), 0 );
    for ( size_t i = 0; i < N.size(); i++ ) {
        if ( N[i] > 0 ) {
            d_submeshes[i].setProcs( N[i] );
            cost[i] = d_maxCostRank;
        }
    }
    setCost0( N, d_submeshes, cost );
    // Add processors until we have used all processors
    int N0 = min( N ) == 0 ? 1 : 0;
    while ( N_proc - sum( N ) - N0 > 0 ) {
        int i = findMax( cost );
        if ( N[i] == 0 ) {
            int k    = -1;
            double c = 0;
            for ( size_t j = 0; j < N.size(); j++ ) {
                if ( N[j] == 0 && c < d_submeshes[i].d_maxCostRank ) {
                    k = j;
                    c = d_submeshes[i].d_maxCostRank;
                }
            }
            N[k] = 1;
            d_submeshes[i].setProcs( 1 );
            cost[i] = d_submeshes[i].d_maxCostRank;
            setCost0( N, d_submeshes, cost );
            N0 = min( N ) == 0 ? 1 : 0;
        } else {
            N[i]++;
            d_submeshes[i].addRank( N[i] - 1 );
            cost[i] = d_submeshes[i].d_maxCostRank;
        }
    }
}
void loadBalanceSimulator::addRank( int rank )
{
    // Add the rank locally
    if ( d_ranks.size() >= d_max_procs )
        return;
    PROFILE_SCOPED( timer, "addRank", 2 );
    for ( int r : d_ranks )
        AMP_ASSERT( r != rank );
    d_ranks.push_back( rank );
    if ( d_submeshes.size() <= 1 || d_ranks.size() == 1 || d_method == 0 ) {
        for ( auto &mesh : d_submeshes )
            mesh.addRank( rank );
    } else if ( ( d_method == 1 || d_method == 2 ) && ( d_ranks.size() > d_submeshes.size() ) ) {
        // Initialize the cost
        std::vector<int> N( d_submeshes.size(), 0 );
        std::vector<double> cost( d_submeshes.size(), 0 );
        for ( size_t i = 0; i < N.size(); i++ ) {
            N[i] = d_submeshes[i].getRanks().size();
            if ( N[i] > 0 )
                cost[i] = d_maxCostRank;
        }
        setCost0( N, d_submeshes, cost );
        // Add the processor to the highest cost mesh
        int i = findMax( cost );
        if ( N[i] == 0 ) {
            d_submeshes[i].setRanks( { rank } );
        } else {
            d_submeshes[i].addRank( rank );
        }
    } else {
        setRanks( d_ranks );
    }
    if ( d_submeshes.empty() )
        d_maxCostRank = d_cost / d_ranks.size();
    else
        d_maxCostRank = max( getRankCost() );
}
void loadBalanceSimulator::setRanks( const std::vector<int> &ranks_in )
{
    PROFILE_SCOPED( timer, "setRanks", 1 );
    // Set the processors for this mesh
    d_ranks = ranks_in;
    AMP_ASSERT( !d_ranks.empty() );
    if ( d_ranks.size() > d_max_procs )
        d_ranks.resize( d_max_procs );
    int N_proc = d_ranks.size();
    int N_mesh = d_submeshes.size();
    // Choose the load balance method
    if ( N_mesh <= 1 || N_proc == 1 || d_method == 0 ) {
        // Everybody is on the same communicator
        for ( auto &mesh : d_submeshes )
            mesh.setRanks( d_ranks );
    } else if ( allEqual( d_submeshes ) ) {
        // Special case where all children are the same size
        if ( N_proc <= N_mesh ) {
            for ( int i = 0; i < N_mesh; i++ ) {
                int j = ( i * N_proc ) / N_mesh;
                d_submeshes[i].setRanks( { d_ranks[j] } );
            }
        } else {
            for ( int i = 0, k = 0, Np = N_proc; i < N_mesh; i++ ) {
                int Ng = N_mesh - i;
                int N  = Np / Ng;
                std::vector<int> ranks2;
                for ( int k2 = 0; k2 < N; k2++, k++ )
                    ranks2.push_back( d_ranks[k] );
                Np -= N;
                d_submeshes[i].setRanks( ranks2 );
            }
        }
    } else if ( d_method == 1 ) {
        // We want to split the meshes onto independent processors
        if ( N_proc == N_mesh ) {
            // Special case where the # of meshes == # of ranks
            for ( int i = 0; i < N_mesh; i++ )
                d_submeshes[i].setRanks( { d_ranks[i] } );
        } else if ( N_proc > N_mesh ) {
            // We have more processors than meshes
            std::vector<int> N( N_mesh, 0 );
            // Perform the initial load balance using 80% of the processors
            int NP = 0.8 * N_proc;
            NP     = std::max( { NP, N_proc - 10 * N_mesh, N_mesh } );
            for ( size_t i = 0; i < N.size(); i++ )
                N[i] = std::max<int>( 1, NP * d_submeshes[i].d_cost / d_cost );
            // Load balance
            loadBalance( N_proc, N );
            AMP_ASSERT( min( N ) > 0 );
            AMP_ASSERT( sum( N ) == N_proc );
            // Set the ranks
            for ( int i = 0, k = 0; i < N_mesh; i++ ) {
                std::vector<int> ranks2;
                for ( int k2 = 0; k2 < N[i]; k2++, k++ )
                    ranks2.push_back( d_ranks[k] );
                d_submeshes[i].setRanks( ranks2 );
            }
        } else {
            // The number of meshes is <= the number of processors
            std::vector<double> cost( N_mesh, 0 );
            for ( int i = 0; i < N_mesh; i++ )
                cost[i] = d_submeshes[i].d_cost;
            auto groups = divideGroups( N_proc, cost );
            for ( int i = 0; i < N_mesh; i++ )
                d_submeshes[i].setRanks( { d_ranks[groups[i]] } );
        }
    } else if ( d_method == 2 ) {
        // Try to achieve the best load balance while splitting the meshes if possible
        if ( N_proc < N_mesh ) {
            // The number of meshes is <= the number of processors
            std::vector<double> cost( N_mesh, 0 );
            for ( int i = 0; i < N_mesh; i++ )
                cost[i] = d_submeshes[i].d_cost;
            auto groups = divideGroups( N_proc, cost );
            for ( int i = 0; i < N_mesh; i++ )
                d_submeshes[i].setRanks( { d_ranks[groups[i]] } );
        } else {
            // Perform the initial load balance using 80% of the processors
            std::vector<int> N( N_mesh, 0 );
            int NP = 0.8 * ( N_proc - 1 );
            NP     = std::max( { NP, N_proc - 10 * N_mesh } );
            for ( size_t i = 0; i < N.size(); i++ )
                N[i] = NP * d_submeshes[i].d_cost / d_cost;
            // Load balance
            loadBalance( N_proc, N );
            if ( min( N ) == 0 )
                AMP_ASSERT( sum( N ) == N_proc - 1 );
            else
                AMP_ASSERT( sum( N ) == N_proc );
            // Set the ranks
            for ( int i = 0, k = 0; i < N_mesh; i++ ) {
                if ( N[i] == 0 ) {
                    d_submeshes[i].setRanks( { d_ranks.back() } );
                } else {
                    std::vector<int> ranks2;
                    for ( int k2 = 0; k2 < N[i]; k2++, k++ )
                        ranks2.push_back( d_ranks[k] );
                    d_submeshes[i].setRanks( ranks2 );
                }
            }
        }
    } else {
        AMP_ERROR( "Unknown method" );
    }
    if ( d_submeshes.empty() )
        d_maxCostRank = d_cost / N_proc;
    else
        d_maxCostRank = max( getRankCost() );
}


} // namespace AMP::Mesh
