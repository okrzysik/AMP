#include "AMP/IO/RestartManager.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/MeshParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/vectors/VectorBuilder.h"


#include <complex>


void record( AMP::UnitTest &ut, bool pass, const std::string &msg )
{
    if ( pass )
        ut.passes( msg );
    else
        ut.failure( msg );
}


void testRestartManager( AMP::UnitTest &ut, const std::string &input_file )
{
    // Load input file
    auto db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database = db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( AMP_COMM_WORLD );

    // Create the meshes from the input database
    auto mesh = AMP::Mesh::MeshFactory::create( params );

    // Create a vector
    auto DOFs = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, 1, 1, true );
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "vec" );
    auto vec = AMP::LinearAlgebra::createVector( DOFs, var, true );
    vec->setToScalar( 100 );

    // Create the restart manager and register data
    AMP::IO::RestartManager writer;
    std::complex<double> complex( 15.3, 4.2 );
    writer.registerData( "int", 13 );
    writer.registerData( "double", 15.3 );
    writer.registerData( "complex", complex );
    writer.registerData( "string", "test" );
    writer.registerData( "inputs", db );
    writer.registerData( "mesh", mesh );
    // writer.registerData( "vec", vec );

    // Write the restart data
    writer.write( "testRestartData" );

    // Read and check the restart data
    AMP::IO::RestartManager reader( "testRestartData" );
    auto int2 = *( reader.getData<int>( "int" ) );
    auto dbl  = *( reader.getData<double>( "double" ) );
    auto cmp  = *( reader.getData<std::complex<double>>( "complex" ) );
    auto str  = *( reader.getData<std::string>( "string" ) );
    bool pass = int2 == 13 && dbl == 15.3 && cmp == complex && str == "test";
    record( ut, pass, "Basic check" );
    auto db2 = reader.getData<AMP::Database>( "inputs" );
    record( ut, db2 != nullptr, "Database loaded" );
    if ( db2 )
        record( ut, *db == *db2, "Database match" );
    auto mesh2 = reader.getData<AMP::Mesh::Mesh>( "mesh" );
    record( ut, mesh2 != nullptr, "Mesh loaded" );
    // auto vec2  = reader.getData<AMP::LinearAlgebra::Vector>( fid, "vec" );
    // record( ut, Vec2 != nullptr, "Load Vector" );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    if ( argc != 2 )
        std::cerr << "test_RestartManager <input>\n";
    testRestartManager( ut, argv[1] );

    int N_failed = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    AMP::AMPManager::shutdown();
    return N_failed;
}
