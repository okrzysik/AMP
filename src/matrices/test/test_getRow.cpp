#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/GetRowHelper.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "ProfilerApp.h"

// Create a vector
std::shared_ptr<AMP::Discretization::DOFManager> createDOFs(
    std::shared_ptr<const AMP::Mesh::Mesh> mesh, const AMP::Mesh::GeomType &type, bool multivec )
{
    auto dof1 = AMP::Discretization::simpleDOFManager::create( mesh, type, 1, 1, true );
    if ( !multivec )
        return dof1;
    auto dof2 = AMP::Discretization::simpleDOFManager::create( mesh, type, 1, 3, true );
    std::vector<std::shared_ptr<AMP::Discretization::DOFManager>> dofs = { dof1, dof2 };
    return std::make_shared<AMP::Discretization::multiDOFManager>( AMP_COMM_WORLD, dofs );
}


// Struct to store test times
struct TestTimes {
    double DOFs          = 0;
    double defaultNNZ    = 0;
    double defaultGetRow = 0;
    double nativeNNZ     = 0;
    double nativeGetRow  = 0;
    void print( const std::string_view prefix )
    {
        if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getRank() != 0 )
            return;
        std::cout << prefix << "Time to create DOFManager: " << DOFs << std::endl;
        std::cout << prefix << "Time to get NNZ (default): " << defaultNNZ << std::endl;
        std::cout << prefix << "Time to get rows (default): " << defaultGetRow << std::endl;
        std::cout << prefix << "Time to get NNZ (native): " << nativeNNZ << std::endl;
        std::cout << prefix << "Time to get rows (native): " << nativeGetRow << std::endl;
    }
};


// Test GetRowHelper
void testGetRowHelper( std::shared_ptr<const AMP::Mesh::Mesh> mesh,
                       const AMP::Mesh::GeomType &leftType,
                       bool leftMultiVec,
                       const AMP::Mesh::GeomType &rightType,
                       bool rightMultiVec,
                       AMP::UnitTest &ut )
{
    TestTimes times;
    AMP::AMP_MPI comm( AMP_COMM_WORLD );

    bool pass = true;

    // Create a DOF managers
    auto start_time = AMP::Utilities::time();
    auto left       = createDOFs( mesh, leftType, leftMultiVec );
    auto right      = createDOFs( mesh, rightType, rightMultiVec );
    comm.barrier();
    times.DOFs = AMP::Utilities::time() - start_time;

    // range of rows to query on this rank
    size_t first_row = left->beginDOF();
    size_t last_row  = left->endDOF();
    size_t num_rows  = last_row - first_row;

    // range of columns for determining local vs remote
    size_t first_col = right->beginDOF();
    size_t last_col  = right->endDOF();

    // create the get row functions from MatrixBuilder
    // default for non-native matrices, and pair using
    // GetRowHelper for native matrices
    // auto leftDOF       = left.get();
    // auto rightDOF      = right.get();
    auto defaultGetRow = [left, right]( size_t row ) {
        auto id = left->getElementID( row );
        return right->getRowDOFs( id );
    };
    auto rowHelper       = std::make_shared<AMP::LinearAlgebra::GetRowHelper>( left, right );
    auto nativeGetRowNNZ = [rowHelper]( size_t row, int &nlocal, int &nremote ) {
        rowHelper->NNZ( row, nlocal, nremote );
    };
    auto nativeGetRow = [rowHelper]( size_t row, size_t *clocal, size_t *cremote ) {
        rowHelper->getRow( row, clocal, cremote );
    };

    // Test default implementation for non-zeros per row
    std::vector<int> NNZ_local( num_rows, 0 );
    std::vector<int> NNZ_remote( num_rows, 0 );
    size_t default_num_nnz_local = 0, default_num_nnz_remote = 0;
    start_time = AMP::Utilities::time();
    for ( size_t row = first_row, i = 0; row < last_row; i++, row++ ) {
        auto cols = defaultGetRow( row );
        for ( auto &&col : cols ) {
            if ( first_col <= col && col < last_col ) {
                NNZ_local[i]++;
                ++default_num_nnz_local;
            } else {
                NNZ_remote[i]++;
                ++default_num_nnz_remote;
            }
        }
    }
    comm.barrier();
    times.defaultNNZ = AMP::Utilities::time() - start_time;

    // Test default implementation for column indices
    std::vector<size_t> local_cols( default_num_nnz_local, 0 ),
        remote_cols( default_num_nnz_remote, 0 );
    size_t local_pos = 0, remote_pos = 0;
    start_time = AMP::Utilities::time();
    for ( size_t row = first_row; row < last_row; row++ ) {
        auto cols = defaultGetRow( row );
        for ( auto &&col : cols ) {
            if ( first_col <= col && col < last_col ) {
                local_cols[local_pos] = col;
                ++local_pos;
            } else {
                remote_cols[local_pos] = col;
                ++remote_pos;
            }
        }
    }
    comm.barrier();
    times.defaultGetRow = AMP::Utilities::time() - start_time;

    // clear out backing vectors and retest with native get row routines
    std::fill( NNZ_local.begin(), NNZ_local.end(), 0 );
    std::fill( NNZ_remote.begin(), NNZ_remote.end(), 0 );
    std::fill( local_cols.begin(), local_cols.end(), 0 );
    std::fill( remote_cols.begin(), remote_cols.end(), 0 );

    // Test native implementation for non-zeros per row
    size_t native_num_nnz_local = 0, native_num_nnz_remote = 0;
    start_time = AMP::Utilities::time();
    for ( size_t row = first_row; row < last_row; row++ ) {
        int nl = 0, nr = 0;
        nativeGetRowNNZ( row, nl, nr );
        native_num_nnz_local += nl;
        native_num_nnz_remote += nr;
    }
    comm.barrier();
    times.nativeNNZ = AMP::Utilities::time() - start_time;

    if ( native_num_nnz_local != default_num_nnz_local ||
         native_num_nnz_remote != default_num_nnz_remote ) {
        pass = false;
        ut.failure( "Default and native NNZ counts differ" );
    }

    // Test default implementation for column indices
    local_pos  = 0;
    remote_pos = 0;
    start_time = AMP::Utilities::time();
    for ( size_t row = first_row, i = 0; row < last_row; i++, row++ ) {
        size_t *clocal  = &local_cols[local_pos];
        size_t *cremote = &remote_cols[remote_pos];
        nativeGetRow( row, clocal, cremote );
        local_pos += NNZ_local[i];
        remote_pos += NNZ_remote[i];
    }
    comm.barrier();
    times.nativeGetRow = AMP::Utilities::time() - start_time;

    if ( pass ) {
        ut.passes( "getRowHelper.getRow" );
    }

    times.print( "" );
}


// Run the tests
void testGetRows( const std::string &input_file, AMP::UnitTest &ut )
{
    AMP::AMP_MPI comm( AMP_COMM_WORLD );

    // Load the input database
    auto db = AMP::Database::parseInputFile( input_file );
    db->print( AMP::plog );
    comm.barrier();

    // Create the Mesh
    auto t1 = AMP::Utilities::time();
    AMP_INSIST( db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto mesh = AMP::Mesh::MeshFactory::create( meshParams );
    comm.barrier();
    auto t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to create mesh: " << t2 - t1 << std::endl;

    // Test GetRowHelper
    auto vertex = AMP::Mesh::GeomType::Vertex;
    testGetRowHelper( mesh, vertex, false, vertex, false, ut );
}


int main( int argc, char **argv )
{

    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;
    PROFILE_ENABLE( 2 );

    for ( int i = 1; i < argc; i++ )
        testGetRows( argv[i], ut );

    ut.report();
    PROFILE_SAVE( "test_getRow" );
    int num_failed = ut.NumFailGlobal();
    ut.reset();
    AMP::AMPManager::shutdown();
    return num_failed;
}
