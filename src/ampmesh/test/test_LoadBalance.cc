// This program simulates the load balance with a given input file on a given number of processors

#include "ampmesh/Mesh.h"
#include "utils/AMPManager.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"

#include <iostream>
#include <iomanip>
#include <math.h>


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
    double t0 = AMP::AMP_MPI::time();
    AMP::Mesh::Mesh::simulated_mesh_struct mesh = AMP::Mesh::Mesh::simulateBuildMesh( params, comm_ranks );
    double t1 = AMP::AMP_MPI::time();

    // Print the results of the load balance
    std::vector<size_t> N_elements(N_procs,0);
    countElements( mesh, N_elements );
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
    std::cout << std::endl;
    std::cout << "min = " << mesh.min() << std::endl;
    std::cout << "max = " << mesh.max() << std::endl;
    std::cout << "avg = " << mesh.avg() << std::endl;
    std::cout << "time = " << t1-t0 << std::endl;

    AMP::AMPManager::shutdown();
    return 0;
}
