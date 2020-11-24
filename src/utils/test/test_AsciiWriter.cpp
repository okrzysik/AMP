#include <sstream>
#include <string>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "ProfilerApp.h"

#include "AMP/utils/Writer.h"

#ifdef USE_AMP_MESH
#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#endif
#ifdef USE_AMP_DISCRETIZATION
#include "AMP/discretization/simpleDOF_Manager.h"
#endif
#ifdef USE_AMP_VECTORS
#include "AMP/vectors/VectorBuilder.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/MatrixBuilder.h"
#endif


// Function to build a vector using a mesh
#if defined( USE_AMP_MESH ) && defined( USE_AMP_VECTORS )
template<int SIZE_X, int SIZE_Y, int SIZE_Z>
AMP::LinearAlgebra::Vector::shared_ptr createVector( AMP::LinearAlgebra::Variable::shared_ptr var,
                                                     AMP::AMP_MPI comm )
{
    std::vector<int> size( 3 );
    size[0] = SIZE_X;
    size[1] = SIZE_Y;
    size[2] = SIZE_Z;
    std::vector<double> range( 6, 0.0 );
    range[1] = 1.0;
    range[3] = 1.0;
    range[5] = 1.0;
    // Create a generic MeshParameters object
    auto database = std::make_shared<AMP::Database>( "Mesh" );
    database->putScalar<int>( "dim", 3 );
    database->putScalar<std::string>( "MeshName", "mesh1" );
    database->putScalar<std::string>( "Generator", "cube" );
    database->putVector<int>( "Size", size );
    database->putVector<double>( "Range", range );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( comm );
    // Create an AMP mesh
    auto mesh = AMP::Mesh::BoxMesh::generate( params );
    // Create the DOF Manager
    auto DOF = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    // Create the vector
    return AMP::LinearAlgebra::createVector( DOF, var, true );
}
#endif


void test_AsciiWriter( AMP::UnitTest *ut )
{

    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::AMP_MPI selfComm( AMP_COMM_SELF );

    // Create the ascii writer
    AMP::Utilities::Writer::shared_ptr writer = AMP::Utilities::Writer::buildWriter( "ASCII" );

// Create and register a vector
#ifdef USE_AMP_VECTORS
    std::string rankString = AMP::Utilities::intToString( globalComm.getRank() + 1, 1 );
    auto var1              = std::make_shared<AMP::LinearAlgebra::Variable>( "vec_global" );
    auto var2              = std::make_shared<AMP::LinearAlgebra::Variable>( "vec_" + rankString );
#ifdef USE_AMP_MESH
    auto vec1 = createVector<2, 3, 4>( var1, globalComm );
    auto vec2 = createVector<3, 2, 1>( var2, selfComm );
#else
    auto vec1 = AMP::LinearAlgebra::createSimpleVector<double>( 20, var1, globalComm );
    auto vec2 = AMP::LinearAlgebra::createSimpleVector<double>( 50, var2, selfComm );
#endif
    writer->registerVector( vec1 );
    writer->registerVector( vec2 );
#endif

// Create and register a matrix
#ifdef USE_AMP_MATRICES
    bool test_matrix = true;
#if !defined( USE_EXT_PETSC ) && !defined( USE_EXT_TRILINOS )
    if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getSize() > 1 ) {
        ut->expected_failure( "No parallel matrix to test" );
        test_matrix = false;
    }
#endif
    if ( test_matrix ) {
        auto mat1 = AMP::LinearAlgebra::createMatrix( vec1, vec1 );
        auto mat2 = AMP::LinearAlgebra::createMatrix( vec2, vec2 );
        writer->registerMatrix( mat1 );
        writer->registerMatrix( mat2 );
        mat1->setScalar( 1.0 );
        mat2->setScalar( globalComm.getRank() + 1 );
    }
#endif

    // Write the output file
    writer->setDecomposition( 1 );
    writer->writeFile( "test_AsciiWriter", 0 );

    ut->passes( "test ran to completion" );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    test_AsciiWriter( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
