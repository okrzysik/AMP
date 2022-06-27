#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/matrices/MatrixBuilder.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/libmesh/ReadTestMesh.h"
#include "AMP/mesh/libmesh/initializeLibMesh.h"
#include "AMP/mesh/libmesh/libmeshMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

// libmesh headers
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"

#include <iostream>
#include <string>


void myTest( AMP::UnitTest *ut, std::string mesh_file )
{
    std::string log_file = "output_testMatMultiply";
    AMP::logOnlyNodeZero( log_file );

    auto mesh_file_db = AMP::Database::parseInputFile( mesh_file );

    // Create a libmesh mesh
    AMP::AMP_MPI comm( AMP_COMM_SELF );
    auto libmeshInit  = std::make_shared<AMP::Mesh::initializeLibMesh>( comm );
    uint32_t mesh_dim = 3;
    libMesh::Parallel::Communicator libMeshComm( comm.getCommunicator() );
    auto myMesh = std::make_shared<libMesh::Mesh>( libMeshComm, mesh_dim );
    AMP::readTestMesh( mesh_file, myMesh );
    libMesh::MeshCommunication().broadcast( *( myMesh.get() ) );
    myMesh->prepare_for_use( false );

    // Create the AMP mesh
    auto myMeshAdapter = std::make_shared<AMP::Mesh::libmeshMesh>( myMesh, "myMesh" );

    // Create the DOF manager
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        myMeshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );

    // Create the vectors
    auto myVar = std::make_shared<AMP::LinearAlgebra::Variable>( "myVar" );
    auto vec1  = AMP::LinearAlgebra::createVector( DOFs, myVar );
    vec1->setToScalar( 1.0 );

    // Create the matrix
    auto mat1 = AMP::LinearAlgebra::createMatrix( vec1, vec1 );
    mat1->zero();
    mat1->setDiagonal( vec1 );
    auto mat2 = mat1->cloneMatrix();
    mat2->zero();
    mat2->setDiagonal( vec1 );

    auto mat3 = AMP::LinearAlgebra::Matrix::matMultiply( mat1, mat2 );

    std::vector<size_t> cols1;
    std::vector<double> vals1;
    mat1->getRowByGlobalID( 0, cols1, vals1 );

    std::vector<size_t> cols2;
    std::vector<double> vals2;
    mat2->getRowByGlobalID( 0, cols2, vals2 );

    std::vector<size_t> cols3;
    std::vector<double> vals3;
    mat3->getRowByGlobalID( 0, cols3, vals3 );

    ut->passes( "testMatMultiply" );

    // Free all data relying on mesh
    DOFs.reset();
    vec1.reset();
    mat1.reset();
    mat2.reset();
    mat3.reset();

    // Free the mesh in the proper order
    myMeshAdapter.reset();
    myMesh.reset();
    libmeshInit.reset();
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    myTest( &ut, "mpcMesh-1" );

    ut.report();

    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
