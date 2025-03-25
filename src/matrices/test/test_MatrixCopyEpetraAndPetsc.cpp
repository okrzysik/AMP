#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


void test_MatricCopyEpetraAndPetsc( AMP::UnitTest &ut )
{
    using DOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>>;

#if defined( AMP_USE_PETSC )
    auto NativePetscFactoryDOF1 = std::make_shared<DOF1>( "NativePetscMatrix" );
#endif

#if defined( AMP_USE_TRILINOS )
    auto EpetraFactory = std::make_shared<DOF1>( "ManagedEpetraMatrix" );

    #if defined( AMP_USE_PETSC )
    test_matrix_loop( ut, NativePetscFactoryDOF1, EpetraFactory );
    test_matrix_loop( ut, EpetraFactory, NativePetscFactoryDOF1 );
    // Test the ManagedPetscMatrix -- TODO
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

    test_MatricCopyEpetraAndPetsc( ut );

    ut.report();
    PROFILE_SAVE( "test_MatrixCopyPetscAndCSR" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
