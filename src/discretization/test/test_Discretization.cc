#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/AMP_MPI.h"
#include "test_Discretization.h"
#include "../../ampmesh/test/meshGenerators.h"



// Main function
int main ( int argc , char ** argv )
{
    AMP::AMPManagerProperties startup_properties;
    //startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    // Run the tests
    testSimpleDOFManager<LibMeshCubeGenerator<5> >( &ut );
    testSimpleDOFManager<ExodusReaderGenerator>( &ut );
    testSimpleDOFManager<MultiMeshGenerator>( &ut );
    testMultiDOFManager<LibMeshCubeGenerator<5> >( &ut );
    testMultiDOFManager<MultiMeshGenerator>( &ut );
    testSubsetDOFManager<MultiMeshGenerator,false>( &ut );
    testSubsetDOFManager<MultiMeshGenerator,true>( &ut );

    // Print the results and return
    ut.report ();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
