#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


void test_MatricCopyPetscAndCSR( AMP::UnitTest &ut )
{
    using DOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>>;
    using DOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>>;

    auto CSRFactoryDOF1 = std::make_shared<DOF1>( "CSRMatrix" );
    auto CSRFactoryDOF3 = std::make_shared<DOF3>( "CSRMatrix" );
#if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using libmeshDOF3          = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>>;
    auto CSRLibmeshFactoryDOF3 = std::make_shared<libmeshDOF3>( "CSRMatrix" );
#endif

#if defined( AMP_USE_PETSC )
    auto NativePetscFactoryDOF1 = std::make_shared<DOF1>( "NativePetscMatrix" );
    auto NativePetscFactoryDOF3 = std::make_shared<DOF3>( "NativePetscMatrix" );

    test_matrix_loop( ut, CSRFactoryDOF1, NativePetscFactoryDOF1 );
    test_matrix_loop( ut, CSRFactoryDOF3, NativePetscFactoryDOF3 );
    test_matrix_loop( ut, NativePetscFactoryDOF1, CSRFactoryDOF1 );
    test_matrix_loop( ut, NativePetscFactoryDOF3, CSRFactoryDOF3 );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    auto NativePetscLibmeshFactoryDOF3 = std::make_shared<libmeshDOF3>( "NativePetscMatrix" );
    test_matrix_loop( ut, CSRLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3 );
    test_matrix_loop( ut, NativePetscLibmeshFactoryDOF3, CSRLibmeshFactoryDOF3 );
    #endif
#endif
}


int main( int argc, char **argv )
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    test_MatricCopyPetscAndCSR( ut );

    ut.report();
    PROFILE_SAVE( "test_MatrixCopyPetscAndCSR" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
