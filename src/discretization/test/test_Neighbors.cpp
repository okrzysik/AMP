#include "../../mesh/test/meshGenerators.h"

#include "AMP/IO/PIO.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"

#include <map>

// Main function
int main( int argc, char **argv )
{
    AMP::AMPManagerProperties startup_properties;
    // startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    auto mesh_generator = AMP::unit_test::AMPCubeGenerator<10>();
    mesh_generator.build_mesh();
    auto mesh        = mesh_generator.getMesh();
    const auto &comm = mesh->getComm();
    if ( comm.getSize() > 1 )
        AMP::pout << "Mesh created on multiple ranks" << std::endl;

    // Create a simple DOF manager
    auto vDOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, false );

    std::map<size_t, size_t> vneighbors;

    auto it = mesh->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    for ( ; it != it.end(); ++it ) {
        const auto row = vDOF->getRowDOFs( *it );
        vneighbors[row.size()]++;
    }

    const auto neighborhood_sizes = ( vneighbors.size() == 4 );
    if ( comm.allReduce( neighborhood_sizes ) )
        ut.passes( "Number of types of neighbors for vertex DOFs passes" );
    else {
        ut.failure( "Number of types of neighbors for vertex DOFs fails" );
    }

    // Print the results and return
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
