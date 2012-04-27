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
    FACTORY factory;
    factory.initMesh();
    InstantiateMatrix<FACTORY>::run_test ( ut );
    VerifyGetSetValuesMatrix<FACTORY>::run_test ( ut );
    VerifyAXPYMatrix<FACTORY>::run_test ( ut );
    VerifyScaleMatrix<FACTORY>::run_test ( ut );
    VerifyMultMatrix<FACTORY>::run_test ( ut );
    VerifyGetLeftRightVector<FACTORY>::run_test ( ut );
    VerifyExtractDiagonal<FACTORY>::run_test ( ut );
    factory.endMesh();
}


template <typename FACTORY>
void  test_petsc_matrix_loop ( AMP::UnitTest *ut )
{
}


int main ( int argc , char **argv )
{

    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    test_matrix_loop<SimpleMatrixFactory> ( &ut );

    test_matrix_loop<DOFMatrixTestFactory<3,3,ExodusReaderGenerator<> > > ( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}


