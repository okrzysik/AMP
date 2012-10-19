#include "loadBalance.h"
#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"

#include <iostream>
#include <iomanip>


// inline function to convert two size_t to double, divide, round the result, and return a size_t integer
static inline size_t divide_double(size_t x, size_t y) { 
    return (size_t) floor( (((double)x)/((double)y)) + 0.5 ); 
}


namespace AMP {
namespace Mesh {


/************************************************************
* Constructors                                              *
************************************************************/
LoadBalance::LoadBalance( )
{
    d_N_elements = 0;
    d_decomp = 0;
    d_min = 0;
    d_max = 0;
    cache_valid = false;
}
LoadBalance::LoadBalance( boost::shared_ptr<MeshParameters> params, const std::vector<int> &ranks, size_t N_elements )
{
    // Get required values from the parameters
    AMP_ASSERT(ranks.size()>0);
    boost::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    AMP_INSIST(database->keyExists("MeshName"),"MeshName must exist in input database");
    d_name = database->getString("MeshName");
    d_type = database->getString("MeshType");
    // Simulate the load process
    d_decomp = 0;
    cache_valid = false;
    if ( d_type == std::string("Multimesh") ) {
        LoadBalance rhs = MultiMesh::simulateBuildMesh( params, ranks );
        d_name = rhs.d_name;
        d_type = rhs.d_type;
        d_N_elements = rhs.d_N_elements;
        d_params = rhs.d_params;
        d_ranks = rhs.d_ranks;
        d_submeshes = rhs.d_submeshes;
        d_decomp = rhs.d_decomp;
        cache_valid = rhs.cache_valid;
        if ( d_ranks.size()==1 ) {
            d_min = d_N_elements;
            d_max = d_N_elements;
            cache_valid = true;
        }
    } else {
        d_params = params;
        d_ranks = ranks;
        d_submeshes.resize(0);
        if ( N_elements==0 )
            d_N_elements = Mesh::estimateMeshSize( params );
        else 
            d_N_elements = N_elements;
        d_min = divide_double(d_N_elements,d_ranks.size());
        d_max = divide_double(d_N_elements,d_ranks.size());
        cache_valid = true;
    }
}
LoadBalance::LoadBalance( boost::shared_ptr<MeshParameters> params, const std::vector<int>& ranks, 
    const std::vector<LoadBalance>& submeshes, int decomp )
{
    boost::shared_ptr<AMP::Database> database = params->getDatabase();
    AMP_ASSERT(database!=NULL);
    AMP_INSIST(database->keyExists("MeshType"),"MeshType must exist in input database");
    AMP_INSIST(database->keyExists("MeshName"),"MeshName must exist in input database");
    d_name = database->getString("MeshName");
    d_type = database->getString("MeshType");
    d_params = params;
    d_ranks = ranks;
    d_submeshes = submeshes;
    d_decomp = decomp;
    cache_valid = false;
    d_N_elements = 0;
    for (size_t i=0; i<d_submeshes.size(); i++)
        d_N_elements += d_submeshes[i].d_N_elements;
    if ( d_ranks.size()==1 ) {
        d_min = d_N_elements;
        d_max = d_N_elements;
        cache_valid = true;
    }
}
LoadBalance::LoadBalance( const LoadBalance& rhs ):
    d_name(rhs.d_name),
    d_type(rhs.d_type),
    d_ranks(rhs.d_ranks),
    d_submeshes(rhs.d_submeshes)
{

    d_N_elements = rhs.d_N_elements;
    d_params = rhs.d_params;
    d_decomp = rhs.d_decomp;
    cache_valid = rhs.cache_valid;
    d_min = rhs.d_min;
    d_max = rhs.d_max;
}


/************************************************************
* Function to add a processor                               *
************************************************************/
void LoadBalance::addProc( int rank )
{
    // Add the rank to the local list
    cache_valid = false;
    d_ranks.push_back(rank);
    if ( d_submeshes.empty() ) 
        return;
    // Choose which submesh to add the processor to
    if ( d_type == std::string("Multimesh") ) {
        MultiMesh::addProcSimulation( *this, d_submeshes, rank, d_decomp );
    } else {
        AMP_ERROR("Not finished");
    }
}


/************************************************************
* Function to return the min, max, and avg # of elements    *
************************************************************/
size_t LoadBalance::min() 
{
    if ( !cache_valid )
        updateCache();
    return d_min;
}
size_t LoadBalance::max() 
{
   if ( !cache_valid )
        updateCache();
    return d_max;
}
size_t LoadBalance::avg() 
{
    return divide_double(d_N_elements,d_ranks.size());
}


/************************************************************
* Function to print the load balance                        *
************************************************************/
void LoadBalance::print() 
{
    int N_procs = 0;
    for (size_t i=0; i<d_ranks.size(); i++)
        N_procs = std::max(N_procs,d_ranks[i]+1);
    std::vector<size_t> N_elements(N_procs,0);
    countElements( *this, N_elements );
    std::cout << "Rank, N_elements:" << std::endl;
    int N_line = 16;
    for (int i=0; i<(N_procs+N_line-1)/N_line; i++) {
        for (int j=i*N_line; j<std::min((i+1)*N_line,N_procs); j++)
            std::cout << std::setw(8) << j;
        std::cout << std::endl;
        for (int j=i*N_line; j<std::min((i+1)*N_line,N_procs); j++)
            std::cout << std::setw(8) << N_elements[j];
        std::cout << std::endl << std::endl;
    }
}


/************************************************************
* Misc. functions                                           *
************************************************************/
void LoadBalance::changeRanks( const std::vector<int> &ranks )
{
    if ( d_submeshes.empty() ) {
        d_ranks = ranks;
        updateCache();
    } else {
        AMP_INSIST(ranks.size()==d_ranks.size(),"Cannot change ranks to different size");
        for (size_t i=0; i<ranks.size(); i++)
            d_ranks[i] = ranks[i];
    }
}
void LoadBalance::countElements( const LoadBalance &mesh, std::vector<size_t> &N_elements )
{
    if ( mesh.d_submeshes.empty() ) {
        for (size_t i=0; i<mesh.d_ranks.size(); i++)
            N_elements[mesh.d_ranks[i]] += mesh.d_N_elements/mesh.d_ranks.size();
    } else {
        for (size_t i=0; i<mesh.d_submeshes.size(); i++)
            countElements( mesh.d_submeshes[i], N_elements );
    }
}
void LoadBalance::updateCache() 
{
    if ( d_submeshes.empty() ) {
        d_min = divide_double(d_N_elements,d_ranks.size());
        d_max = divide_double(d_N_elements,d_ranks.size());
        return;
    }
    if ( d_decomp==0 ) {
        // General case
        int N_procs = 0;
        for (size_t i=0; i<d_ranks.size(); i++)
            N_procs = std::max(N_procs,d_ranks[i]+1);
        std::vector<size_t> N_elements(N_procs,0);
        countElements( *this, N_elements );
        d_min = d_N_elements;
        d_max = 0;
        for (size_t i=0; i<d_ranks.size(); i++) {
            d_min = std::min(d_min,N_elements[d_ranks[i]]);
            d_max = std::max(d_max,N_elements[d_ranks[i]]);
        }
    } else if ( d_decomp==1 ) {
        // Special case where no two submeshes share a processor
        d_min = d_N_elements;
        d_max = 0;
        for (size_t i=0; i<d_submeshes.size(); i++) {
            if ( d_submeshes[i].cache_valid ) {
                if ( d_submeshes[i].d_min < d_min )
                    d_min = d_submeshes[i].d_min;
                if ( d_submeshes[i].d_max > d_max )
                    d_max = d_submeshes[i].d_max;
            } else {
                d_min = std::min(d_min,d_submeshes[i].min());
                d_max = std::max(d_max,d_submeshes[i].max());
            }
        }
    } else {
        AMP_ERROR("Unkown value for d_decomp");
    }
    cache_valid = true;
}



}
}

