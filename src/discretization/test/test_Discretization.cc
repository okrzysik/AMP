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
    testSimpleDOFManager<AMPMeshCubeGenerator<10> >( &ut );
    #ifdef USE_LIBMESH
        testSimpleDOFManager<LibMeshCubeGenerator<5> >( &ut );
        testSimpleDOFManager<ExodusReaderGenerator<1> >( &ut );
        testSimpleDOFManager<ExodusReaderGenerator<3> >( &ut );
    #endif
    testSimpleDOFManager<MultiMeshGenerator>( &ut );
    testMultiDOFManager<AMPMeshCubeGenerator<10> >( &ut );
    #ifdef USE_LIBMESH
        testMultiDOFManager<LibMeshCubeGenerator<5> >( &ut );
    #endif
    testMultiDOFManager<MultiMeshGenerator>( &ut );
    testSubsetDOFManager<MultiMeshGenerator,false>( &ut );
    testSubsetDOFManager<MultiMeshGenerator,true>( &ut );
    #ifdef USE_LIBMESH
        testSubsetDOFManager<ExodusReaderGenerator<3>,false>( &ut );
    #else
        testSubsetDOFManager<AMPMeshCubeGenerator<10>,false>( &ut );
    #endif

    // Print the results and return
    ut.report ();
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
