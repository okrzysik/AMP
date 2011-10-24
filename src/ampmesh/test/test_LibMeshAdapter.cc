#include "string.h"
#include "test_MeshAdapterLoop.h"
#include "utils/UnitTest.h"
#include "utils/AMPManager.h"



int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    Bug_623::run_test( &ut );

    if ( ut.size() == 1 ) 
        MeshAdapterLoop <ThreeElementLGenerator>::test_mesh ( &ut );
    else
        ut.expected_failure("ThreeElementLGenerator is meant for serial tests");
    MeshAdapterLoop <LibMeshCubeGenerator<5> >::test_mesh ( &ut );
    MeshAdapterLoop <ExodusReaderGenerator>::test_mesh ( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;

}
