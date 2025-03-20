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
    // Test the CSRMatrix
    test_matrix_loop<CSRFactoryDOF1, CSRFactoryDOF1>( ut );
    test_matrix_loop<CSRFactoryDOF3, CSRFactoryDOF3>( ut );

#if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    test_matrix_loop<CSRLibmeshFactoryDOF3, CSRLibmeshFactoryDOF3>( ut );
#endif

#if defined( AMP_USE_PETSC )
    using NativePetscFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 3>;
    using NativePetscFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 3>;
#endif

#if defined( AMP_USE_TRILINOS )
    using EpetraFactory = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>;
    test_matrix_loop<EpetraFactory, EpetraFactory>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using EpetraLibmeshFactory = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 1>;
    test_matrix_loop<EpetraLibmeshFactory, EpetraLibmeshFactory>( ut );
    #endif
#endif

    // Test the DenseSerialMatrix (only valid for serial meshes)
    // Note: we need to be careful, the memory requirement can be very large
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() == 1 ) {
        using DenseSerialFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 2>;
        using DenseSerialFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 2>;
        test_matrix_loop<DenseSerialFactoryDOF1, DenseSerialFactoryDOF1>( ut );
        test_matrix_loop<DenseSerialFactoryDOF3, DenseSerialFactoryDOF3>( ut );
    }

#if defined( AMP_USE_PETSC )

    test_matrix_loop<NativePetscFactoryDOF1, NativePetscFactoryDOF1>( ut );
    test_matrix_loop<NativePetscFactoryDOF3, NativePetscFactoryDOF3>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using NativePetscLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 3>;
    test_matrix_loop<NativePetscLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3>( ut );
    #endif
#endif

    ut.report();
    PROFILE_SAVE( "test_MatrixCopy" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
