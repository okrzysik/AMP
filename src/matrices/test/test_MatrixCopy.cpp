#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


void testCopy( AMP::UnitTest &ut )
{
    using DOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>>;
    using DOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>>;

    auto CSRFactoryDOF1 = std::make_shared<DOF1>( "CSRMatrix" );
    auto CSRFactoryDOF3 = std::make_shared<DOF3>( "CSRMatrix" );
#if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using libmeshDOF3          = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>>;
    auto CSRLibmeshFactoryDOF3 = std::make_shared<libmeshDOF3>( "CSRMatrix" );
#endif
    // Test the CSRMatrix
    test_matrix_loop( ut, CSRFactoryDOF1, CSRFactoryDOF1 );
    test_matrix_loop( ut, CSRFactoryDOF3, CSRFactoryDOF3 );

#if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    test_matrix_loop( ut, CSRLibmeshFactoryDOF3, CSRLibmeshFactoryDOF3 );
#endif

#if defined( AMP_USE_PETSC )
    auto NativePetscFactoryDOF1 = std::make_shared<DOF1>( "NativePetscMatrix" );
    auto NativePetscFactoryDOF3 = std::make_shared<DOF3>( "NativePetscMatrix" );
#endif

#if defined( AMP_USE_TRILINOS )
    auto EpetraFactory = std::make_shared<DOF1>( "ManagedEpetraMatrix" );
    test_matrix_loop( ut, EpetraFactory, EpetraFactory );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    auto EpetraLibmeshFactory = std::make_shared<libmeshDOF3>( "ManagedEpetraMatrix" );
    test_matrix_loop( ut, EpetraLibmeshFactory, EpetraLibmeshFactory );
    #endif
#endif

    // Test the DenseSerialMatrix (only valid for serial meshes)
    // Note: we need to be careful, the memory requirement can be very large
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() == 1 ) {
        auto DenseSerialFactoryDOF1 = std::make_shared<DOF1>( "DenseSerialMatrix" );
        auto DenseSerialFactoryDOF3 = std::make_shared<DOF1>( "DenseSerialMatrix" );
        test_matrix_loop( ut, DenseSerialFactoryDOF1, DenseSerialFactoryDOF1 );
        test_matrix_loop( ut, DenseSerialFactoryDOF3, DenseSerialFactoryDOF3 );
    }

#if defined( AMP_USE_PETSC )

    test_matrix_loop( ut, NativePetscFactoryDOF1, NativePetscFactoryDOF1 );
    test_matrix_loop( ut, NativePetscFactoryDOF3, NativePetscFactoryDOF3 );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    auto NativePetscLibmeshFactoryDOF3 = std::make_shared<libmeshDOF3>( "NativePetscMatrix" );
    test_matrix_loop( ut, NativePetscLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3 );
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

    testCopy( ut );

    ut.report();
    PROFILE_SAVE( "test_MatrixCopy" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
