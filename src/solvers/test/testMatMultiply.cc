#include <iostream>
#include <string>

#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/ReadTestMesh.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "matrices/MatrixBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"

// libmesh headers
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_data.h"
#include "libmesh/mesh_generation.h"


void myTest( AMP::UnitTest *ut, std::string mesh_file )
{
    std::string log_file = "output_testMatMultiply";
    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> mesh_file_db( new AMP::InputDatabase( "mesh_file_db" ) );
    AMP::InputManager::getManager()->parseInputFile( mesh_file, mesh_file_db );

    // Create a libmesh mesh
    AMP::AMP_MPI comm( AMP_COMM_SELF );
    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( comm ) );
    const unsigned int mesh_dim = 3;
    AMP::shared_ptr<::Mesh> myMesh( new ::Mesh( mesh_dim ) );
    AMP::readTestMesh( mesh_file, myMesh );
    MeshCommunication().broadcast( *( myMesh.get() ) );
    myMesh->prepare_for_use( false );

    // Create the AMP mesh
    AMP::Mesh::Mesh::shared_ptr myMeshAdapter( new AMP::Mesh::libMesh( myMesh, "myMesh" ) );

    // Create the DOF manager
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( myMeshAdapter, AMP::Mesh::Vertex, 1, 3 );

    // Create the vectors
    AMP::LinearAlgebra::Variable::shared_ptr myVar( new AMP::LinearAlgebra::Variable( "myVar" ) );
    AMP::LinearAlgebra::Vector::shared_ptr vec1 = AMP::LinearAlgebra::createVector( DOFs, myVar );
    vec1->setToScalar( 1.0 );

    // Create the matrix
    AMP::LinearAlgebra::Matrix::shared_ptr mat1 = AMP::LinearAlgebra::createMatrix( vec1, vec1 );
    mat1->zero();
    mat1->setDiagonal( vec1 );
    AMP::LinearAlgebra::Matrix::shared_ptr mat2 = mat1->cloneMatrix();
    mat2->zero();
    mat2->setDiagonal( vec1 );

    AMP::LinearAlgebra::Matrix::shared_ptr mat3 =
        AMP::LinearAlgebra::Matrix::matMultiply( mat1, mat2 );

    std::vector<unsigned int> cols1;
    std::vector<double> vals1;
    mat1->getRowByGlobalID( 0, cols1, vals1 );

    std::vector<unsigned int> cols2;
    std::vector<double> vals2;
    mat2->getRowByGlobalID( 0, cols2, vals2 );

    std::vector<unsigned int> cols3;
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
