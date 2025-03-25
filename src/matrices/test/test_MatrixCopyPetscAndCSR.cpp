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

    using CSRFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 4>;
    using CSRFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 4>;
#if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using CSRLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 4>;
#endif

#if defined( AMP_USE_PETSC )
    using NativePetscFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 3>;
    using NativePetscFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 3>;

    test_matrix_loop<CSRFactoryDOF1, NativePetscFactoryDOF1>( ut );
    test_matrix_loop<CSRFactoryDOF3, NativePetscFactoryDOF3>( ut );
    test_matrix_loop<NativePetscFactoryDOF1, CSRFactoryDOF1>( ut );
    test_matrix_loop<NativePetscFactoryDOF3, CSRFactoryDOF3>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using NativePetscLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 3>;
    test_matrix_loop<CSRLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3>( ut );
    test_matrix_loop<NativePetscLibmeshFactoryDOF3, CSRLibmeshFactoryDOF3>( ut );
    #endif
#endif

    ut.report();
    PROFILE_SAVE( "test_MatrixCopyPetscAndCSR" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
