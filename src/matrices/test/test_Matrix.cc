#include "vectors/trilinos/EpetraVectorEngine.h"
#include "test_Matrix.h"
#include "test_MatrixTests.h"
#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"


using namespace AMP::unit_test;


template <typename FACTORY>
void  test_matrix_loop ( AMP::UnitTest *ut )
{
    InstantiateMatrix<FACTORY>::run_test ( ut );
    VerifyGetSetValuesMatrix<FACTORY>::run_test ( ut );
    VerifyAXPYMatrix<FACTORY>::run_test ( ut );
    VerifyScaleMatrix<FACTORY>::run_test ( ut );
    VerifyMultMatrix<FACTORY>::run_test ( ut );
    VerifyGetLeftRightVector<FACTORY>::run_test ( ut );
    // VerifyExtractDiagonal<FACTORY>::run_test ( ut );
}


template <typename FACTORY>
void  test_petsc_matrix_loop ( AMP::UnitTest *ut )
{
}


int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    SimpleMatrixFactory::initMesh ();
    test_matrix_loop<SimpleMatrixFactory> ( &ut );
    SimpleMatrixFactory::endMesh();

    MeshMatrixTestFactory<3,3,ExodusReaderGenerator>::initMesh();
    test_matrix_loop<MeshMatrixTestFactory<3,3,ExodusReaderGenerator> > ( &ut );
    MeshMatrixTestFactory<3,3,ExodusReaderGenerator>::endMesh();

    // MeshMatrixTestFactory< 1 , 3 , ExodusReaderGenerator>::initMesh();
    // test_matrix_loop<MeshMatrixTestFactory < 1 , 3 , ExodusReaderGenerator> > ( ut );
    // MeshMatrixTestFactory< 1 , 3 , ExodusReaderGenerator>::endMesh();

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


