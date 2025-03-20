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

#if defined( AMP_USE_TRILINOS )
    using EpetraFactory = DOFMatrixTestFactory<1, 1, AMPCubeGenerator<5>, 1>;
    test_matrix_loop<CSRFactoryDOF1, EpetraFactory>( ut );
    test_matrix_loop<EpetraFactory, CSRFactoryDOF1>( ut );

    #if defined( AMP_USE_LIBMESH ) && defined( USE_AMP_DATA )
    using EpetraLibmeshFactory  = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 1>;
    using CSRLibmeshFactoryDOF3 = DOFMatrixTestFactory<3, 3, ExodusReaderGenerator<>, 4>;
    test_matrix_loop<EpetraLibmeshFactory, CSRLibmeshFactoryDOF3>( ut );
    test_matrix_loop<CSRLibmeshFactoryDOF3, EpetraLibmeshFactory>( ut );
    #endif

#endif

    ut.report();
    PROFILE_SAVE( "test_MatrixCopyEpetraAndCSR" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
