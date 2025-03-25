#include "test_Matrix.h"

#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/matrices/testHelpers/MatrixTests.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/vectors/VectorBuilder.h"


using namespace AMP::LinearAlgebra;


void test_MatricCopyEpetraAndCSR( AMP::UnitTest &ut )
{
    using DOF1 = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>>;

    auto CSRFactoryDOF1 = std::make_shared<DOF1>( "CSRMatrix" );

#if defined( AMP_USE_TRILINOS )
    auto EpetraFactory = std::make_shared<DOF1>( "ManagedEpetraMatrix" );
    test_matrix_loop( ut, CSRFactoryDOF1, EpetraFactory );
    test_matrix_loop( ut, EpetraFactory, CSRFactoryDOF1 );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using libmeshFactory       = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>>;
    auto EpetraLibmeshFactory  = std::make_shared<libmeshFactory>( "ManagedEpetraMatrix" );
    auto CSRLibmeshFactoryDOF3 = std::make_shared<libmeshFactory>( "CSRMatrix" );
    test_matrix_loop( ut, EpetraLibmeshFactory, CSRLibmeshFactoryDOF3 );
    test_matrix_loop( ut, CSRLibmeshFactoryDOF3, EpetraLibmeshFactory );
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

    test_MatricCopyEpetraAndCSR( ut );

    ut.report();
    PROFILE_SAVE( "test_MatrixCopyPetscAndCSR" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
