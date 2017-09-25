#include "test_Matrix.h"
#include "ProfilerApp.h"
#include "test_MatrixTests.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "matrices/MatrixBuilder.h"
#include "vectors/SimpleVector.h"


using namespace AMP::LinearAlgebra;


template <typename FACTORY>
void test_matrix_loop( AMP::UnitTest &ut )
{
    std::string name = FACTORY::name();
    PROFILE_START( name );
    FACTORY factory;
    factory.initMesh();
    InstantiateMatrix<FACTORY>::run_test( &ut );
    VerifyGetSetValuesMatrix<FACTORY>::run_test( &ut );
    VerifyAXPYMatrix<FACTORY>::run_test( &ut );
    VerifyScaleMatrix<FACTORY>::run_test( &ut );
    VerifyGetLeftRightVector<FACTORY>::run_test( &ut );
    VerifyExtractDiagonal<FACTORY>::run_test( &ut );
    VerifyMultMatrix<FACTORY>::run_test( &ut );
    VerifyMatMultMatrix<FACTORY>::run_test( &ut );
    VerifyAddElementNode<FACTORY>::run_test( &ut );
    factory.endMesh();
    PROFILE_STOP( name );
}


template <typename FACTORY>
void test_petsc_matrix_loop( AMP::UnitTest &ut )
{
    NULL_USE( ut );
}


void testBasics( AMP::UnitTest &ut, const std::string &type )
{
    // Test creating a non-square matrix and ensure it is the proper size
    auto left  = AMP::LinearAlgebra::SimpleVector<double>::create( 5, "left" );
    auto right = AMP::LinearAlgebra::SimpleVector<double>::create( 10, "right" );
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
#ifdef USE_EXT_TRILINOS
    types.push_back( "ManagedEpetraMatrix" );
#endif
#if defined( USE_EXT_TRILINOS ) && defined( USE_EXT_PETSC )
    types.push_back( "ManagedPetscMatrix" );
#endif
    types.push_back( "auto" );
    for ( auto type : types )
        testBasics( ut, type );

// Test the ManagedPetscMatrix (requires both petsc and trilinos)
#if defined( USE_EXT_PETSC ) && defined( USE_EXT_TRILINOS )
    test_matrix_loop<DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>>( ut );
    test_matrix_loop<DOFMatrixTestFactory<3, 3, AMPCubeGenerator<5>, 1>>( ut );
#ifdef USE_EXT_LIBMESH
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
    AMP::AMPManager::shutdown();
    return num_failed;
}
