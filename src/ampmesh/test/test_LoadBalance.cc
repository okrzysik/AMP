// This program simulates the load balance with a given input file on a given number of processors

#include "ampmesh/Mesh.h"
#include "utils/AMPManager.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"

#include <iostream>
#include <iomanip>


// Function to recursively cound the number of elements per processor
void countElements( const AMP::Mesh::Mesh::simulated_mesh_struct &mesh, std::vector<size_t> &N_elements )
{
    if ( mesh.submeshes.empty() ) {
        for (size_t i=0; i<mesh.ranks.size(); i++)
            N_elements[mesh.ranks[i]] += mesh.N_elements/mesh.ranks.size();
    } else {
        for (size_t i=0; i<mesh.submeshes.size(); i++)
            countElements( mesh.submeshes[i], N_elements );
    }
}


// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManager::startup(argc,argv);

    // Load the input file and the desired number of processors
    if (argc != 3) {
        std::cout << "Error calling test_LoadBalancer, format should be:" << std::endl;
        std::cout << "   ./test_LoadBalancer  N_procs  input_file" << std::endl;
        return -1;
    }
    int N_procs = std::atoi(argv[1]);
    std::string input_file = argv[2];

    // Simulate loading the mesh
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file, input_db );
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    std::vector<int> comm_ranks(N_procs);
    for (int i=0; i<N_procs; i++)
        comm_ranks[i] = i;
    AMP::Mesh::Mesh::simulated_mesh_struct mesh = AMP::Mesh::Mesh::simulateBuildMesh( params, comm_ranks );

    // Check the results of the load balance
    std::vector<size_t> N_elements(N_procs,0);
    countElements( mesh, N_elements );
    std::cout << "Rank   N_elements" << std::endl;
    for (int i=0; i<N_procs; i++)
        std::cout << std::setw(4) << i << std::setw(10) << N_elements[i] << std::endl;

    AMP::AMPManager::shutdown();
    return 0;
}
