#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


int main( int argc, char **argv )
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    PROFILE_ENABLE();

#if defined( AMP_USE_PETSC )
    using NativePetscFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 3>;
#endif

#if defined( AMP_USE_TRILINOS )
    using EpetraFactory = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>;

    #if defined( AMP_USE_PETSC )
    test_matrix_loop<NativePetscFactoryDOF1, EpetraFactory>( ut );
    test_matrix_loop<EpetraFactory, NativePetscFactoryDOF1>( ut );
    // Test the ManagedPetscMatrix -- TODO
    #endif
#endif

    ut.report();
    PROFILE_SAVE( "test_MatrixCopy" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
