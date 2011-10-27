#include "string.h"

#include "test_MeshAdapterLoop.h"

#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"



int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    if ( ut.size() == 1 ) 
        MeshAdapterVectorLoop <ThreeElementLGenerator>::test_mesh ( &ut );
    else
        ut.expected_failure("ThreeElementLGenerator is meant for serial tests");
    MeshAdapterVectorLoop <LibMeshCubeGenerator<5> >::test_mesh ( &ut );
    MeshAdapterVectorLoop <ExodusReaderGenerator>::test_mesh ( &ut );

    ut.report ();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
