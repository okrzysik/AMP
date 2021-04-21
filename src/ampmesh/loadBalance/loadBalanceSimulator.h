#ifndef included_AMP_loadBalanceSimulator
#define included_AMP_loadBalanceSimulator

#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/utils/Database.h"


namespace AMP::Mesh {


/**
 * \class loadBalanceSimulator
 * \brief Class used to contain simulated mesh load
 * \details  This structure provides info that can be used to simulate loading
 *   a mesh, and checking the resulting load balance
 */
class loadBalanceSimulator
{
public:
    /**
     * \brief    Default constructor
     * \details  This will simulate creating a new load balance
     * \param params  Input parameters for the mesh that will be used
     * \param ranks   List of processor ranks that will be used
     * \param N_elements    Optional argument specifying the number of elements on the mesh
     *                      (0: Get the number of elements through a call to
     * Mesh::estimateMeshSize())
     */
    loadBalanceSimulator( std::shared_ptr<AMP::Database> params );

    /**
     * \brief    Default constructor
     * \details  This will simulate creating a new load balance
     * \param cost      The cost of the mesh (# of elements)
     * \param maxProc   The maximum number of processors (0: no limit)
     * \param name      The mesh name
     */
    loadBalanceSimulator( double cost, int maxProc = 0, const std::string &name = "" );

    /**
     * \brief    Default constructor
     * \details  This will simulate creating a new load balance
     * \param meshes    The list of meshes
     * \param method    The load balance method
     * \param name      The mesh name
     */
    loadBalanceSimulator( const std::vector<loadBalanceSimulator> &meshes,
                          int method              = 2,
                          const std::string &name = "" );

    //! Empty constructor
    loadBalanceSimulator();

    //! Function to get the current ranks in the load balance
    size_t getMethod() const { return d_method; }

    //! Function to get the ranks for the mesh
    const std::vector<int> &getRanks() const { return d_ranks; }

    //! Function to get the ranks for the ith submesh
    inline const std::vector<int> &getRanks( int i ) const { return d_submeshes[i].d_ranks; }

    //! Function to change the ranks
    void setProcs( int N_ranks );

    //! Function to change the ranks
    void setRanks( const std::vector<int> &ranks );

    //! Return the submeshes
    inline const auto &getSubmeshes() const { return d_submeshes; }

    //! Function to get the total cost
    double getCost() const { return d_cost; }

    //! Get the cost per rank ( # of elements * relative weight )
    std::vector<double> getRankCost() const;

    /**
     * \brief    Print the mesh hierarchy
     * \details  This function will print the load balance and mesh hierarchy
     * \param detail    The details on what to print (bit array)
     *                  Bit 0: print the load balance by rank
     *                  Bit 1: print the number of procs per mesh
     * \param indent    Number of spaces to indent the printing
     */
    void print( uint8_t detail = 3, uint8_t indent = 0 );


private:
    // Internal data
    std::string d_name;
    double d_cost;
    double d_maxCostRank;
    size_t d_max_procs;
    int d_method;
    std::vector<int> d_ranks;
    std::vector<loadBalanceSimulator> d_submeshes;

    void addRank( int rank );

    // Function to add the rank cost
    void addRankCost( std::vector<double> &cost ) const;
    void loadBalance( int, std::vector<int> &N );
};


} // namespace AMP::Mesh

#endif
