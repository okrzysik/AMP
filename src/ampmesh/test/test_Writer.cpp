#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/MultiMesh.h"
#include "AMP/ampmesh/structured/BoxMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#ifdef USE_AMP_VECTORS
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#endif
#ifdef USE_AMP_MATRICES
#include "AMP/matrices/MatrixBuilder.h"
#endif

#include "ProfilerApp.h"

#include <sstream>
#include <string>


// Get the surface element type given the volume type
inline AMP::Mesh::GeomType getSurfaceType( AMP::Mesh::GeomType volume )
{
    AMP_ASSERT( volume != AMP::Mesh::GeomType::null );
    if ( volume == AMP::Mesh::GeomType::Vertex )
        return AMP::Mesh::GeomType::Vertex;
    return static_cast<AMP::Mesh::GeomType>( static_cast<int>( volume ) - 1 );
}


// Calculate the volume of each element
AMP::LinearAlgebra::Vector::shared_ptr calcVolume( AMP::Mesh::Mesh::shared_ptr mesh )
{
    auto DOF =
        AMP::Discretization::simpleDOFManager::create( mesh, mesh->getGeomType(), 0, 1, false );
    auto var = std::make_shared<AMP::LinearAlgebra::Variable>( "volume" );
    auto vec = AMP::LinearAlgebra::createVector( DOF, var, true );
    vec->zero();
    std::vector<size_t> dofs;
    for ( const auto &elem : mesh->getIterator( mesh->getGeomType(), 0 ) ) {
        double volume = elem.volume();
        DOF->getDOFs( elem.globalID(), dofs );
        AMP_ASSERT( dofs.size() == 1 );
        vec->addLocalValuesByGlobalID( 1, &dofs[0], &volume );
    }
    return vec;
}


// Print the mesh names
void printMeshNames( AMP::Mesh::Mesh::shared_ptr mesh, const std::string &prefix = "" )
{
    std::cout << prefix << mesh->getName() << std::endl;
    auto multimesh = std::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( mesh );
    if ( multimesh ) {
        for ( auto mesh2 : multimesh->getMeshes() )
            printMeshNames( mesh2, prefix + "   " );
    }
}
void printMeshNames( const std::string &filename )
{
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    auto input_db = AMP::Database::parseInputFile( filename );
    input_db->print( AMP::plog );
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );
    auto mesh = AMP::Mesh::Mesh::buildMesh( params );
    if ( globalComm.getRank() == 0 ) {
        std::cout << "Mesh names (rank 0):\n";
        printMeshNames( mesh, "   " );
    }
}


// Function to build a vector using a mesh
#if defined( USE_AMP_MESH ) && defined( USE_AMP_VECTORS )
template<int SIZE_X, int SIZE_Y, int SIZE_Z>
AMP::LinearAlgebra::Vector::shared_ptr createVector( AMP::LinearAlgebra::Variable::shared_ptr var,
                                                     AMP::AMP_MPI comm )
{
    // Create an AMP mesh
    auto size     = { SIZE_X, SIZE_Y, SIZE_Z };
    auto range    = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
    auto database = AMP::Database::create(
        "dim", 3, "MeshName", "mesh", "Generator", "cube", "Size", size, "Range", range );
    auto params = std::make_shared<AMP::Mesh::MeshParameters>( std::move( database ) );
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


/***************************************************************
 * Test the writer with a simple vector                         *
 ***************************************************************/
void testWriterVector( AMP::UnitTest &ut, const std::string &writerName )
{
#ifdef USE_AMP_VECTORS

    // Create the writer and get it's properties
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    auto writer     = AMP::Utilities::Writer::buildWriter( writerName, comm );
    auto properties = writer->getProperties();
    if ( !properties.registerVector ) {
        ut.expected_failure( writerName + " does not support registering an independent vector" );
        return;
    }

    // Create and register the vectors
    AMP::ArraySize block( 0, 0, comm.getRank() );
    auto var1 = std::make_shared<AMP::LinearAlgebra::Variable>( "vec1" );
    auto var2 = std::make_shared<AMP::LinearAlgebra::Variable>( "vec2" );
    auto var3 = std::make_shared<AMP::LinearAlgebra::Variable>( "vec3" );
    auto vec1 = AMP::LinearAlgebra::createArrayVector<double>( { 5, 4, 3 }, block, comm, var1 );
    auto vec2 = AMP::LinearAlgebra::createSimpleVector<double>( 20, var2, comm );
    auto vec3 = AMP::LinearAlgebra::createSimpleVector<double>( 50, var3, AMP_COMM_SELF );
    vec1->setRandomValues();
    vec2->setRandomValues();
    vec3->setRandomValues();
    writer->registerVector( vec1, "vec1" );
    writer->registerVector( vec2, "vec2" );
    writer->registerVector( vec3, "vec3" );

    // Write the file
    auto rankStr         = std::to_string( comm.getRank() + 1 );
    std::string filename = "output_test_Writer/vector-" + writerName + "-" + rankStr + "proc";
    writer->writeFile( filename, 0, 0.0 );
    if ( AMP::Utilities::fileExists( filename + "_0." + properties.extension ) )
        ut.passes( writerName + " registered independent vector" );
    else
        ut.failure( writerName + " registered independent vector" );

#endif
}


/***************************************************************
 * Test the writer with a simple matrix                         *
 ***************************************************************/
void testWriterMatrix( AMP::UnitTest &ut, const std::string &writerName )
{
#ifdef USE_AMP_MATRICES

    // Create the writer and get it's properties
    AMP::AMP_MPI comm( AMP_COMM_WORLD );
    auto writer     = AMP::Utilities::Writer::buildWriter( writerName, comm );
    auto properties = writer->getProperties(); // Function to build a vector using a mesh
    if ( !properties.registerMatrix ) {
        ut.expected_failure( writerName + " does not support registering a matrix" );
        return;
    }
#if !defined( USE_EXT_PETSC ) && !defined( USE_EXT_TRILINOS )
    ut.expected_failure( writerName + "  - no parallel matrix to test" );
    return;
#endif

    // Create and register a matrix
    auto rankStr = std::to_string( comm.getRank() + 1 );
    auto var1    = std::make_shared<AMP::LinearAlgebra::Variable>( "vec_global" );
    auto var2    = std::make_shared<AMP::LinearAlgebra::Variable>( "vec_" + rankStr );
    auto vec1    = createVector<2, 3, 4>( var1, comm );
    auto vec2    = createVector<3, 2, 1>( var2, AMP_COMM_SELF );
    auto mat1    = AMP::LinearAlgebra::createMatrix( vec1, vec1 );
    auto mat2    = AMP::LinearAlgebra::createMatrix( vec2, vec2 );
    writer->registerVector( vec1, "vec1" );
    writer->registerVector( vec2, "vec2" );
    writer->registerMatrix( mat1 );
    writer->registerMatrix( mat2 );
    mat1->setScalar( 1.0 );
    mat2->setScalar( comm.getRank() + 1 );

    // Write the file
    std::string filename = "output_test_Writer/matrix-" + writerName + "-" + rankStr + "proc";
    writer->setDecomposition( 1 );
    writer->writeFile( filename, 0, 0.0 );
    if ( AMP::Utilities::fileExists( filename + "_0." + properties.extension ) )
        ut.passes( writerName + " registered matrix" );
    else
        ut.failure( writerName + " registered matrix" );
#endif
}


/***************************************************************
 * Test the writer with an input file                           *
 ***************************************************************/
void testWriterMesh( AMP::UnitTest &ut,
                     const std::string &writerName,
                     const std::string &input_file )
{
    // Create the writer and get it's properties
    auto writer     = AMP::Utilities::Writer::buildWriter( writerName, AMP_COMM_WORLD );
    auto properties = writer->getProperties();

    // Check that we have valid work to do
    if ( !properties.registerMesh ) {
        ut.expected_failure( writerName + " does not support registering a mesh" );
        return;
    }

    // Start the timer
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    globalComm.barrier();
    double t1 = AMP::AMP_MPI::time();

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    PROFILE_START( "Load Mesh" );
    auto mesh        = AMP::Mesh::Mesh::buildMesh( params );
    auto pointType   = AMP::Mesh::GeomType::Vertex;
    auto volumeType  = mesh->getGeomType();
    auto surfaceType = getSurfaceType( volumeType );
    globalComm.barrier();
    PROFILE_STOP( "Load Mesh" );
    double t2 = AMP::AMP_MPI::time();

    // Create a surface mesh
    auto submesh = mesh->Subset( mesh->getSurfaceIterator( surfaceType, 1 ) );

#ifdef USE_AMP_VECTORS
    // Create a simple DOFManager
    auto DOFparams  = std::make_shared<AMP::Discretization::DOFManagerParameters>( mesh );
    auto DOF_scalar = AMP::Discretization::simpleDOFManager::create( mesh, pointType, 1, 1, true );
    auto DOF_volume = AMP::Discretization::simpleDOFManager::create( mesh, volumeType, 1, 1, true );
    auto DOF_vector = AMP::Discretization::simpleDOFManager::create( mesh, pointType, 1, 3, true );
    auto DOF_gauss  = AMP::Discretization::simpleDOFManager::create( mesh, volumeType, 1, 8, true );

    // Create the vectors
    auto rank_var     = std::make_shared<AMP::LinearAlgebra::Variable>( "rank" );
    auto position_var = std::make_shared<AMP::LinearAlgebra::Variable>( "position" );
    auto gp_var       = std::make_shared<AMP::LinearAlgebra::Variable>( "gp_var" );
    auto id_var       = std::make_shared<AMP::LinearAlgebra::Variable>( "ids" );
    auto meshID_var   = std::make_shared<AMP::LinearAlgebra::Variable>( "MeshID" );
    auto block_var    = std::make_shared<AMP::LinearAlgebra::Variable>( "block" );
    auto rank_vec     = AMP::LinearAlgebra::createVector( DOF_scalar, rank_var, true );
    auto position     = AMP::LinearAlgebra::createVector( DOF_vector, position_var, true );
    auto gauss_pt     = AMP::LinearAlgebra::createVector( DOF_gauss, gp_var, true );
    auto meshID_vec   = AMP::LinearAlgebra::createVector( DOF_scalar, meshID_var, true );
    auto block_vec    = AMP::LinearAlgebra::createVector( DOF_volume, block_var, true );
    gauss_pt->setToScalar( 100 );
    globalComm.barrier();
#endif
    double t3 = AMP::AMP_MPI::time();

// Create a view of a vector
#ifdef USE_AMP_VECTORS
    AMP::Mesh::Mesh::shared_ptr clad;
    AMP::LinearAlgebra::Vector::shared_ptr z_surface;
    AMP::LinearAlgebra::Vector::shared_ptr cladPosition;
    if ( submesh ) {
        AMP::LinearAlgebra::VS_MeshIterator meshSelector( submesh->getIterator( pointType, 1 ),
                                                          submesh->getComm() );
        AMP::LinearAlgebra::VS_Stride zSelector( 2, 3 );
        auto vec_meshSubset = position->select( meshSelector, "mesh subset" );
        AMP_ASSERT( vec_meshSubset );
        z_surface = vec_meshSubset->select( zSelector, "z surface" );
        AMP_ASSERT( z_surface );
        clad = mesh->Subset( "clad" );
        if ( clad ) {
            clad->setName( "clad" );
            AMP::LinearAlgebra::VS_Mesh cladMeshSelector( clad );
            cladPosition = position->select( cladMeshSelector, "cladPosition" );
        }
    }
#endif

    // Register the data
    int level = 1; // How much detail do we want to register
    writer->registerMesh( mesh, level );
    if ( submesh )
        writer->registerMesh( submesh, level );
#ifdef USE_AMP_VECTORS
    if ( properties.registerVectorWithMesh ) {
        writer->registerVector( meshID_vec, mesh, pointType, "MeshID" );
        writer->registerVector( block_vec, mesh, volumeType, "BlockID" );
        writer->registerVector( rank_vec, mesh, pointType, "rank" );
        writer->registerVector( position, mesh, pointType, "position" );
        writer->registerVector( gauss_pt, mesh, volumeType, "gauss_pnt" );
        if ( submesh )
            writer->registerVector( z_surface, submesh, pointType, "z_surface" );
        // Register a vector over the clad
        if ( clad )
            writer->registerVector( cladPosition, clad, pointType, "clad_position" );
    }
#endif
    globalComm.barrier();
    double t4 = AMP::AMP_MPI::time();

    // For each submesh, store the mesh id, surface ids, and volume
    auto meshIDs = mesh->getBaseMeshIDs();
    for ( size_t i = 0; i < meshIDs.size(); i++ ) {
        auto mesh2 = mesh->Subset( meshIDs[i] );
        if ( mesh2 ) {
            auto volume = calcVolume( mesh2 );
            AMP::LinearAlgebra::VS_Mesh meshSelector( mesh2 );
            auto meshID_vec2 = meshID_vec->select( meshSelector, "mesh subset" );
            meshID_vec2->setToScalar( i + 1 );
            writer->registerMesh( mesh2, level );
            if ( properties.registerVectorWithMesh )
                writer->registerVector( volume, mesh2, mesh2->getGeomType(), "volume" );
            // Get the surface
            auto surfaceMesh = mesh2->Subset( mesh2->getSurfaceIterator( surfaceType, 1 ) );
            if ( surfaceMesh ) {
                auto DOF_surface = AMP::Discretization::simpleDOFManager::create(
                    surfaceMesh, surfaceType, 0, 1, true );
                auto id_vec = AMP::LinearAlgebra::createVector( DOF_surface, id_var, true );
                if ( properties.registerVectorWithMesh )
                    writer->registerVector( id_vec, surfaceMesh, surfaceType, "surface_ids" );
                id_vec->setToScalar( -1 );
                std::vector<size_t> dofs;
                for ( auto &id : surfaceMesh->getBoundaryIDs() ) {
                    double val = double( id );
                    for ( auto elem : surfaceMesh->getBoundaryIDIterator( surfaceType, id, 0 ) ) {
                        DOF_surface->getDOFs( elem.globalID(), dofs );
                        AMP_ASSERT( dofs.size() == 1 );
                        id_vec->setValuesByGlobalID( 1, &dofs[0], &val );
                    }
                }
            }
        }
    }

// Initialize the data
#ifdef USE_AMP_VECTORS
    rank_vec->setToScalar( globalComm.getRank() );
    rank_vec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    std::vector<size_t> dofs;
    for ( auto elem : DOF_vector->getIterator() ) {
        DOF_vector->getDOFs( elem.globalID(), dofs );
        auto pos = elem.coord();
        position->setValuesByGlobalID( dofs.size(), dofs.data(), pos.data() );
    }
    position->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    block_vec->setToScalar( -1 );
    for ( auto &id : mesh->getBlockIDs() ) {
        double val = double( id );
        for ( auto elem : mesh->getBlockIDIterator( volumeType, id, 0 ) ) {
            DOF_volume->getDOFs( elem.globalID(), dofs );
            block_vec->setValuesByGlobalID( 1, &dofs[0], &val );
        }
    }
    block_vec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    globalComm.barrier();
#endif
    double t5 = AMP::AMP_MPI::time();

    // Write a single output file
    std::string fname = "output_test_Writer/" + input_file + "-" + writerName + "-" +
                        std::to_string( globalComm.getSize() );
    if ( globalComm.getSize() <= 20 ) {
        globalComm.barrier();
        writer->setDecomposition( 1 );
        writer->writeFile( fname + "proc_single", 0 );
        globalComm.barrier();
    }
    double t6 = AMP::AMP_MPI::time();

    // Write a separate output file for each rank
    {
        globalComm.barrier();
        writer->setDecomposition( 2 );
        writer->writeFile( fname + "proc_multiple", 0 );
        globalComm.barrier();
    }
    double t7 = AMP::AMP_MPI::time();

    if ( globalComm.getRank() == 0 ) {
        std::cout << writerName << ":" << std::endl;
        std::cout << "  Create meshes: " << t2 - t1 << std::endl;
        std::cout << "  Allocate vectors: " << t3 - t2 << std::endl;
        std::cout << "  Register data: " << t4 - t3 << std::endl;
        std::cout << "  Initialize vectors: " << t5 - t4 << std::endl;
        if ( globalComm.getSize() <= 20 )
            std::cout << "  Write a single file: " << t6 - t5 << std::endl;
        std::cout << "  Write multiple files: " << t7 - t6 << std::endl;
        std::cout << "  Total time: " << t7 - t1 << std::endl;
    }

    ut.passes( writerName + " test ran to completion" );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE();
    AMP::logOnlyNodeZero( "output_test_SiloIO" );

    // const auto writers = { "Silo", "HDF5", "Ascii" };
    const auto writers = { "HDF5" };

    if ( argc == 1 ) {

        // Run basic tests
        for ( auto writer : writers ) {
            testWriterVector( ut, writer );
            testWriterMatrix( ut, writer );
            testWriterMesh( ut, writer, "input_SiloIO-1" );
        }

    } else {

        // Test the provided input files
        for ( int i = 1; i < argc; i++ ) {

            if ( AMP::AMP_MPI( AMP_COMM_WORLD ).getRank() == 0 )
                std::cout << "Testing " << argv[i] << std::endl;

            // Print the mesh names (rank 0)
            printMeshNames( argv[i] );

            // Run the tests
            for ( auto writer : writers )
                testWriterMesh( ut, writer, argv[i] );
        }
    }

    int N_failed = ut.NumFailGlobal();
    ut.report();
    ut.reset();
    PROFILE_SAVE( "test_Silo" );
    AMP::AMPManager::shutdown();
    return N_failed;
}
