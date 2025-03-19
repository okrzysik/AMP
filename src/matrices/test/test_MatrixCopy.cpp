#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


void test_matrix_loop( AMP::UnitTest &ut, std::shared_ptr<MatrixTests> tests )
{
    tests->InstantiateMatrix( &ut );
    tests->VerifyGetSetValuesMatrix( &ut );
    tests->VerifyAXPYMatrix( &ut );
    tests->VerifyCopyMatrix( &ut );
    tests->VerifyScaleMatrix( &ut );
    tests->VerifyGetLeftRightVector( &ut );
    tests->VerifyExtractDiagonal( &ut );
    tests->VerifyMultMatrix( &ut );
    tests->VerifyMatMultMatrix( &ut );
    tests->VerifyAddElementNode( &ut );
}

template<typename FACTORY1, typename FACTORY2>
void test_matrix_loop( AMP::UnitTest &ut )
{
    auto factory          = std::make_shared<FACTORY1>();
    std::string name      = factory->name();
    auto copy_factory     = std::make_shared<FACTORY2>();
    std::string copy_name = copy_factory->name();
    PROFILE2( name + copy_name );
    auto tests = std::make_shared<MatrixTests>( factory, copy_factory );
    test_matrix_loop( ut, tests );
}

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
    test_matrix_loop<CSRFactoryDOF1, EpetraFactory>( ut );
    test_matrix_loop<EpetraFactory, CSRFactoryDOF1>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using EpetraLibmeshFactory = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 1>;
    test_matrix_loop<EpetraLibmeshFactory, EpetraLibmeshFactory>( ut );
    test_matrix_loop<EpetraLibmeshFactory, CSRLibmeshFactoryDOF3>( ut );
    test_matrix_loop<CSRLibmeshFactoryDOF3, EpetraLibmeshFactory>( ut );
    #endif

    #if defined( AMP_USE_PETSC )
    test_matrix_loop<NativePetscFactoryDOF1, EpetraFactory>( ut );
    test_matrix_loop<EpetraFactory, NativePetscFactoryDOF1>( ut );
    // Test the ManagedPetscMatrix -- TODO
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
    test_matrix_loop<CSRFactoryDOF1, NativePetscFactoryDOF1>( ut );
    test_matrix_loop<CSRFactoryDOF3, NativePetscFactoryDOF3>( ut );
    test_matrix_loop<NativePetscFactoryDOF1, CSRFactoryDOF1>( ut );
    test_matrix_loop<NativePetscFactoryDOF3, CSRFactoryDOF3>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using NativePetscLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 3>;
    test_matrix_loop<NativePetscLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3>( ut );
    test_matrix_loop<CSRLibmeshFactoryDOF3, NativePetscLibmeshFactoryDOF3>( ut );
    test_matrix_loop<NativePetscLibmeshFactoryDOF3, CSRLibmeshFactoryDOF3>( ut );
    #endif
#endif

    ut.report();
    PROFILE_SAVE( "test_MatrixCopy" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
