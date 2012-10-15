#include "test_Matrix.h"
#include "test_MatrixTests.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/ProfilerApp.h"


using namespace AMP::unit_test;



template <typename FACTORY>
void  test_matrix_loop ( AMP::UnitTest *ut )
{
    std::string name = "test_matrix_loop: "+FACTORY::name();
    PROFILE_START(name);
    FACTORY factory;
    factory.initMesh();
    #if defined(USE_EXT_PETSC) && defined(USE_EXT_PETSC)
        InstantiateMatrix<FACTORY>::run_test ( ut );
        VerifyGetSetValuesMatrix<FACTORY>::run_test ( ut );
        VerifyAXPYMatrix<FACTORY>::run_test ( ut );
        VerifyScaleMatrix<FACTORY>::run_test ( ut );
        VerifyGetLeftRightVector<FACTORY>::run_test ( ut );
        VerifyExtractDiagonal<FACTORY>::run_test ( ut );
        VerifyMultMatrix<FACTORY>::run_test ( ut );
        VerifyMatMultMatrix<FACTORY>::run_test ( ut );
    #else
        ut->failure("Tests require petsc and trilinos");
    #endif
    factory.endMesh();
    PROFILE_STOP(name);
}


template <typename FACTORY>
void  test_petsc_matrix_loop ( AMP::UnitTest *ut )
{
}


int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    test_matrix_loop<SimpleMatrixFactory> ( &ut );

    #ifdef USE_EXT_LIBMESH
        test_matrix_loop<DOFMatrixTestFactory<3,3,ExodusReaderGenerator<> > > ( &ut );
    #else
        test_matrix_loop<DOFMatrixTestFactory<3,3,AMPCubeGenerator<5> > > ( &ut );
    #endif

    ut.report();
    PROFILE_SAVE("test_Matrix");
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


