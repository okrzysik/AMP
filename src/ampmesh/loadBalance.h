#ifndef included_AMP_LoadBalance
#define included_AMP_LoadBalance

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshParameters.h"


namespace AMP {
namespace Mesh {


/**
 * \class LoadBalance
 * \brief Class used to contain simulated mesh load
 * \details  This structure provides info that can be used to simulate loading 
 *   a mesh, and checking the resulting load balance
 */
class LoadBalance
{
public:

    /**
     * \brief    Default constructor
     * \details  This will simulate creating a new load balance
     * \param params  Input parameters for the mesh that will be used
     * \param ranks   List of processor ranks that will be used
     * \param N_elements    Optional argument specifying the number of elements on the mesh
     *                      (0: Get the number of elements through a call to Mesh::estimateMeshSize())
     */
    LoadBalance( AMP::shared_ptr<MeshParameters> params, const std::vector<int> &ranks, size_t N_elements=0 );

    /**
     * \brief    Advanced constructor
     * \details  This is an advanced constuctor used to create the load balance
     * \param params  Input parameters for the current mesh
     * \param ranks   List of processor ranks that are used
     * \param submeshes  List of sub-meshes
     * \param decomp  Domain decomposition for the submeshes
     *                (0: General, 1: No set of submeshes share a rank, and all ranks are used)
     */
    LoadBalance( AMP::shared_ptr<MeshParameters> params, const std::vector<int>& ranks, const std::vector<LoadBalance>& submeshes, int decomp );

    //! Empty constructor
    LoadBalance();

    //! Copy constructor
    LoadBalance(const LoadBalance&);

    /**
     * \brief    Function to add a processor
     * \details  This function will add a processor with the given rank to the load balance,
     *    adjusting the load balance of any subsequent meshes
     * \param rank  The rank of the processor to add
     */
    bool addProc( int rank );

    //! Function to get the current ranks in the load balance
    size_t getSize() const { return d_N_elements; }

    //! Function to get the total number of elements in the load balance
    const std::vector<int>& getRanks() const { return d_ranks; }

    //! Function to change the ranks
    void changeRanks( const std::vector<int> &ranks );

    //! Function to get the mesh parameters
    AMP::shared_ptr<MeshParameters> getParams() const { return d_params; }

    /**
     * \brief    Print the mesh hierarchy
     * \details  This function will print the load balance and mesh hierarchy
     * \param detail    The details on what to print (bit array)
     *                  Bit 0: print the load balance by rank
     *                  Bit 1: print the number of procs per mesh
     * \param indent    Number of spaces to indent the printing
     */
    void print( unsigned char detail=3, unsigned char indent=0 );

    //! Function to return the minimum number of elements on a processor
    size_t min();

    //! Function to return the maximum number of elements on a processor
    size_t max();

    //! Function to return the average number of elements per processor
    size_t avg();

private:

    // Internal data
    std::string         d_name;
    std::string         d_type;
    size_t              d_N_elements;
    size_t              d_max_ranks;
    std::vector<int>    d_ranks;
    AMP::shared_ptr<MeshParameters> d_params;
    std::vector<LoadBalance>  d_submeshes;

    // Special flag used to identify key decompositions
    char d_decomp;      // 0: General decomposition, 1: All submeshes are on distinct sets of processors

    // Cache some commonly used data
    bool cache_valid;
    size_t d_min, d_max;

    // Function to count the elements on each rank
    void countElements( const LoadBalance &mesh, std::vector<size_t> &N_elements );

    // Function to update cached data
    void updateCache();
};


}
}
#endif

