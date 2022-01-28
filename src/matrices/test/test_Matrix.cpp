#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


template<typename FACTORY>
void test_matrix_loop( AMP::UnitTest &ut )
{
    auto factory     = std::make_shared<FACTORY>();
    std::string name = factory->name();
    PROFILE_START( name );
    MatrixTests tests( factory );
    tests.InstantiateMatrix( &ut );
    tests.VerifyGetSetValuesMatrix( &ut );
    tests.VerifyAXPYMatrix( &ut );
    tests.VerifyScaleMatrix( &ut );
    tests.VerifyGetLeftRightVector( &ut );
    tests.VerifyExtractDiagonal( &ut );
    tests.VerifyMultMatrix( &ut );
    tests.VerifyMatMultMatrix( &ut );
    tests.VerifyAddElementNode( &ut );
    PROFILE_STOP( name );
}


template<typename FACTORY>
void test_petsc_matrix_loop( AMP::UnitTest &ut )
{
    NULL_USE( ut );
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
    std::vector<std::string> types = { "DenseSerialMatrix" };
#ifdef AMP_USE_TRILINOS
    types.emplace_back( "ManagedEpetraMatrix" );
#endif
    types.emplace_back( "auto" );
    for ( auto type : types )
        testBasics( ut, type );

// Test the ManagedPetscMatrix (requires both petsc and trilinos)
#if defined( AMP_USE_PETSC ) && defined( AMP_USE_TRILINOS )
    test_matrix_loop<DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>>( ut );
    test_matrix_loop<DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 1>>( ut );
    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    test_matrix_loop<DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 1>>( ut );
    #endif
#endif

    // Test the DenseSerialMatrix (only valid for serial meshes)
    // Note: we need to be careful, the memory requirement can be very large
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() == 1 ) {
        test_matrix_loop<DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 2>>( ut );
        test_matrix_loop<DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 2>>( ut );
    }

    ut.report();
    PROFILE_SAVE( "test_Matrix" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
