#include "AMP/discretization/DOF_Manager.h"
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


// Get the number of local and remote columns for a given row
std::pair<size_t, size_t> getRowNNZ( const AMP::Discretization::DOFManager *leftDOF,
                                     const AMP::Discretization::DOFManager *rightDOF,
                                     size_t row )
{
    auto elem       = leftDOF->getElement( row );
    auto dofs       = rightDOF->getRowDOFs( elem );
    size_t N_local  = 0;
    size_t N_remote = 0;
    size_t start    = rightDOF->beginDOF();
    size_t end      = rightDOF->endDOF();
    for ( auto dof : dofs ) {
        if ( dof >= start && dof < end )
            N_local++;
        else
            N_remote++;
    }
    return { N_local, N_remote };
}


// Get the local and remote columns for a given row
void getRow( const AMP::Discretization::DOFManager *leftDOF,
             const AMP::Discretization::DOFManager *rightDOF,
             size_t row,
             size_t *local,
             size_t *remote )
{
    auto elem    = leftDOF->getElement( row );
    auto dofs    = rightDOF->getRowDOFs( elem );
    size_t start = rightDOF->beginDOF();
    size_t end   = rightDOF->endDOF();
    size_t i     = 0;
    size_t j     = 0;
    for ( auto dof : dofs ) {
        if ( dof >= start && dof < end )
            local[i++] = dof;
        else
            remote[j++] = dof;
    }
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

    // Create a DOF manager for a nodal vector
    t1               = AMP::Utilities::time();
    auto nodalDofMap = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto leftDOF  = nodalDofMap.get();
    auto rightDOF = nodalDofMap.get();
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to create DOFManager: " << t2 - t1 << std::endl;

    // Test naive implementation for getRowNNZ
    size_t begin = leftDOF->beginDOF();
    size_t end   = leftDOF->endDOF();
    std::vector<size_t> NNZ_local( end - begin, 0 );
    std::vector<size_t> NNZ_remote( end - begin, 0 );
    t1 = AMP::Utilities::time();
    for ( size_t row = begin, i = 0; row < end; i++, row++ )
        std::tie( NNZ_local[i], NNZ_remote[i] ) = getRowNNZ( leftDOF, rightDOF, row );
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to get NNZ (naive): " << t2 - t1 << std::endl;

    // Test naive implementation for getRow
    size_t NNZ[2] = { 0, 0 };
    for ( size_t i = 0; i < NNZ_local.size(); i++ ) {
        NNZ[0] += NNZ_local[i];
        NNZ[1] += NNZ_remote[i];
    }
    auto local     = new size_t[NNZ[0]];
    auto remote    = new size_t[NNZ[1]];
    auto localPtr  = new size_t *[NNZ_local.size()];
    auto remotePtr = new size_t *[NNZ_local.size()];
    for ( size_t i = 0, j1 = 0, j2 = 0; i < NNZ_local.size(); i++ ) {
        localPtr[i]  = &local[j1];
        remotePtr[i] = &remote[j2];
        j1 += NNZ_local[i];
        j2 += NNZ_remote[i];
    }
    t1 = AMP::Utilities::time();
    for ( size_t row = begin, i = 0; row < end; i++, row++ )
        getRow( leftDOF, rightDOF, row, localPtr[i], remotePtr[i] );
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to get rows (naive): " << t2 - t1 << std::endl;

    // Test getRowHelper
    t1 = AMP::Utilities::time();
    AMP::LinearAlgebra::GetRowHelper rowHelper( nodalDofMap, nodalDofMap );
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to construct (GetRowHelper): " << t2 - t1 << std::endl;
    bool pass = true;
    t1        = AMP::Utilities::time();
    for ( size_t row = begin, i = 0; row < end; i++, row++ ) {
        auto [N1, N2] = rowHelper.NNZ( row );
        pass          = pass && N1 == NNZ_local[i] && N2 == NNZ_remote[i];
    }
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to get NNZ (GetRowHelper): " << t2 - t1 << std::endl;
    if ( pass )
        ut.passes( "getRowHelper.NNZ" );
    else
        ut.failure( "getRowHelper.NNZ" );
    pass = true;
    t1   = AMP::Utilities::time();
    size_t d1[1024], d2[1024];
    for ( size_t row = begin, i = 0; row < end; i++, row++ ) {
        rowHelper.getRow( row, d1, d2 );
        for ( size_t j = 0; j < NNZ_local[i]; j++ )
            pass = pass && d1[j] == localPtr[i][j];
        for ( size_t j = 0; j < NNZ_remote[i]; j++ )
            pass = pass && d2[j] == remotePtr[i][j];
    }
    comm.barrier();
    t2 = AMP::Utilities::time();
    if ( comm.getRank() == 0 )
        std::cout << "Time to get row (GetRowHelper): " << t2 - t1 << std::endl;
    if ( pass )
        ut.passes( "getRowHelper.getRow" );
    else
        ut.failure( "getRowHelper.getRow" );


    delete[] local;
    delete[] remote;
    delete[] localPtr;
    delete[] remotePtr;
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
