#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


template<typename FACTORY>
void test_petsc_matrix_loop( [[maybe_unused]] AMP::UnitTest &ut )
{
}

void testBasics( AMP::UnitTest &ut, const std::string &type )
{
    // Test creating a non-square matrix and ensure it is the proper size
    auto left  = AMP::LinearAlgebra::createSimpleVector<double>( 5, "left" );
    auto right = AMP::LinearAlgebra::createSimpleVector<double>( 10, "right" );
    auto mat   = AMP::LinearAlgebra::createMatrix(
        right, left, type, []( size_t row ) { return std::vector<size_t>( 1, row ); } );
    int rows = mat->numGlobalRows();
    int cols = mat->numGlobalColumns();
    if ( rows == 5 && cols == 10 )
        ut.passes( "Created non-square matrix (" + type + ")" );
    else
        ut.failure( "Failed non-square matrix (" + type + ")" );
}


int main( int argc, char **argv )
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    // Test some basic properties
    std::vector<std::string> types = { "DenseSerialMatrix", "CSRMatrix" };
#ifdef AMP_USE_TRILINOS
    types.emplace_back( "ManagedEpetraMatrix" );
#endif
#ifdef AMP_USE_PETSC
    types.emplace_back( "NativePetscMatrix" );
#endif
    types.emplace_back( "auto" );
    for ( auto type : types )
        testBasics( ut, type );

    using CSRFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 4>;
    using CSRFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 4>;
    // Test the CSRMatrix
    test_matrix_loop<CSRFactoryDOF1>( ut );
    test_matrix_loop<CSRFactoryDOF3>( ut );

#if defined( AMP_USE_PETSC )
    using NativePetscFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 3>;
    using NativePetscFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 3>;
#endif

#if defined( AMP_USE_TRILINOS )
    using EpetraFactory = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>;
    test_matrix_loop<EpetraFactory>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using EpetraLibmeshFactory = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 1>;
    test_matrix_loop<EpetraLibmeshFactory>( ut );
    #endif
#endif

    // Test the DenseSerialMatrix (only valid for serial meshes)
    // Note: we need to be careful, the memory requirement can be very large
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() == 1 ) {
        using DenseSerialFactoryDOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 2>;
        using DenseSerialFactoryDOF3 = DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 2>;

        test_matrix_loop<DenseSerialFactoryDOF1>( ut );
        test_matrix_loop<DenseSerialFactoryDOF3>( ut );
    }

#if defined( AMP_USE_PETSC )

    test_matrix_loop<NativePetscFactoryDOF1>( ut );
    test_matrix_loop<NativePetscFactoryDOF3>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using NativePetscLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 3>;
    test_matrix_loop<NativePetscLibmeshFactoryDOF3>( ut );
    #endif
#endif

    ut.report();
    PROFILE_SAVE( "test_Matrix" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
