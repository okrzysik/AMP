#include <sstream>
#include <string>

#include "AMP/ampmesh/Mesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
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

#include "ProfilerApp.h"


AMP::Mesh::GeomType getSurfaceType( AMP::Mesh::GeomType volume )
{
    if ( volume == AMP::Mesh::GeomType::Vertex )
        return AMP::Mesh::GeomType::Vertex;
    else if ( volume == AMP::Mesh::GeomType::Edge )
        return AMP::Mesh::GeomType::Vertex;
    else if ( volume == AMP::Mesh::GeomType::Face )
        return AMP::Mesh::GeomType::Edge;
    else if ( volume == AMP::Mesh::GeomType::Volume )
        return AMP::Mesh::GeomType::Face;
    return AMP::Mesh::GeomType::null;
}


AMP::LinearAlgebra::Vector::shared_ptr calcVolume( AMP::Mesh::Mesh::shared_ptr mesh )
{
    auto DOF =
        AMP::Discretization::simpleDOFManager::create( mesh, mesh->getGeomType(), 0, 1, false );
    auto var = AMP::make_shared<AMP::LinearAlgebra::Variable>( "volume" );
    auto vec = AMP::LinearAlgebra::createVector( DOF, var, true );
    vec->zero();
    std::vector<size_t> dofs;
    for ( const auto &elem : mesh->getIterator( mesh->getGeomType(), 0 ) ) {
        double volume = elem.volume();
        DOF->getDOFs( elem.globalID(), dofs );
        AMP_ASSERT( dofs.size() == 1 );
        vec->addLocalValueByGlobalID( dofs[0], volume );
    }
    return vec;
}


void test_Silo( AMP::UnitTest *ut, std::string input_file )
{

    AMP::PIO::logOnlyNodeZero( "output_test_SiloIO" );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    globalComm.barrier();
    double t1 = AMP::AMP_MPI::time();

    // Read the input file
    auto input_db = AMP::make_shared<AMP::InputDatabase>( "input_db" );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    auto database = input_db->getDatabase( "Mesh" );
    auto params   = AMP::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( globalComm );

    // Create the meshes from the input database
    PROFILE_START( "Load Mesh" );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );
    auto pointType                   = AMP::Mesh::GeomType::Vertex;
    auto volumeType                  = mesh->getGeomType();
    auto surfaceType                 = getSurfaceType( volumeType );
    globalComm.barrier();
    PROFILE_STOP( "Load Mesh" );
    double t2 = AMP::AMP_MPI::time();

    // Create a surface mesh
    auto submesh = mesh->Subset( mesh->getSurfaceIterator( surfaceType, 1 ) );

#ifdef USE_AMP_VECTORS
    // Create a simple DOFManager
    auto DOFparams = AMP::make_shared<AMP::Discretization::DOFManagerParameters>( mesh );
    auto DOF_scalar = AMP::Discretization::simpleDOFManager::create( mesh, pointType, 1, 1, true );
    auto DOF_vector = AMP::Discretization::simpleDOFManager::create( mesh, pointType, 1, 3, true );
    auto DOF_gauss  = AMP::Discretization::simpleDOFManager::create( mesh, volumeType, 1, 8, true );
    auto DOF_surface =
        AMP::Discretization::simpleDOFManager::create( submesh, surfaceType, 0, 1, true );

    // Create the vectors
    auto rank_var     = AMP::make_shared<AMP::LinearAlgebra::Variable>( "rank" );
    auto position_var = AMP::make_shared<AMP::LinearAlgebra::Variable>( "position" );
    auto gp_var       = AMP::make_shared<AMP::LinearAlgebra::Variable>( "gp_var" );
    auto id_var       = AMP::make_shared<AMP::LinearAlgebra::Variable>( "ids" );
    auto rank_vec     = AMP::LinearAlgebra::createVector( DOF_scalar, rank_var, true );
    auto position     = AMP::LinearAlgebra::createVector( DOF_vector, position_var, true );
    auto gauss_pt     = AMP::LinearAlgebra::createVector( DOF_gauss, gp_var, true );
    auto id_vec       = AMP::LinearAlgebra::createVector( DOF_surface, id_var, true );
    gauss_pt->setToScalar( 100 );
    globalComm.barrier();
#endif
    double t3 = AMP::AMP_MPI::time();

// Create a view of a vector
#ifdef USE_AMP_VECTORS
    AMP::Mesh::Mesh::shared_ptr clad;
    AMP::LinearAlgebra::Vector::shared_ptr z_surface;
    AMP::LinearAlgebra::Vector::shared_ptr cladPosition;
    if ( submesh != nullptr ) {
        AMP::LinearAlgebra::VS_MeshIterator meshSelector( submesh->getIterator( pointType, 1 ),
                                                          submesh->getComm() );
        AMP::LinearAlgebra::VS_Stride zSelector( 2, 3 );
        auto vec_meshSubset = position->select( meshSelector, "mesh subset" );
        AMP_ASSERT( vec_meshSubset.get() != nullptr );
        z_surface = vec_meshSubset->select( zSelector, "z surface" );
        AMP_ASSERT( z_surface.get() != nullptr );
        clad = mesh->Subset( "clad" );
        if ( clad.get() != nullptr ) {
            clad->setName( "clad" );
            AMP::LinearAlgebra::VS_Mesh cladMeshSelector( clad );
            cladPosition = position->select( cladMeshSelector, "cladPosition" );
        }
    }
#endif

    // Create the silo writer and register the data
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    int level       = 1; // How much detail do we want to register
    siloWriter->registerMesh( mesh, level );
    if ( submesh != nullptr )
        siloWriter->registerMesh( submesh, level );
#ifdef USE_AMP_VECTORS
    siloWriter->registerVector( rank_vec, mesh, pointType, "rank" );
    siloWriter->registerVector( position, mesh, pointType, "position" );
    siloWriter->registerVector( gauss_pt, mesh, volumeType, "gauss_pnt" );
    if ( submesh != nullptr ) {
        siloWriter->registerVector( z_surface, submesh, pointType, "z_surface" );
        siloWriter->registerVector( id_vec, submesh, surfaceType, "surface_ids" );
    }
    // Register a vector over the clad
    if ( clad.get() != nullptr )
        siloWriter->registerVector( cladPosition, clad, pointType, "clad_position" );
#endif
    globalComm.barrier();
    double t4 = AMP::AMP_MPI::time();

    // For each submesh, store the volume
    auto meshIDs = mesh->getBaseMeshIDs();
    for ( auto id : meshIDs ) {
        auto mesh2 = mesh->Subset( id );
        if ( mesh2 != nullptr ) {
            auto volume = calcVolume( mesh2 );
            siloWriter->registerMesh( mesh2, level );
            siloWriter->registerVector( volume, mesh2, mesh2->getGeomType(), "volume" );
        }
    }

// Initialize the data
#ifdef USE_AMP_VECTORS
    rank_vec->setToScalar( globalComm.getRank() );
    rank_vec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    std::vector<size_t> dofs;
    for ( auto it = DOF_vector->getIterator(); it != it.end(); ++it ) {
        AMP::Mesh::MeshElementID id = it->globalID();
        DOF_vector->getDOFs( id, dofs );
        std::vector<double> pos = it->coord();
        position->setValuesByGlobalID( dofs.size(), &dofs[0], &pos[0] );
    }
    position->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    if ( submesh != nullptr ) {
        id_vec->setToScalar( -1 );
        auto ids = submesh->getBoundaryIDs();
        for ( auto &id : ids ) {
            auto it = submesh->getBoundaryIDIterator( surfaceType, id, 0 );
            for ( size_t j = 0; j < it.size(); j++ ) {
                DOF_surface->getDOFs( it->globalID(), dofs );
                AMP_ASSERT( dofs.size() == 1 );
                id_vec->setValueByGlobalID( dofs[0], id );
                ++it;
            }
        }
    }
    globalComm.barrier();
#endif
    double t5 = AMP::AMP_MPI::time();

    // Write a single output file
    if ( globalComm.getSize() <= 20 ) {
        std::stringstream fname1;
        fname1 << input_file << "_" << globalComm.getSize() << "proc_single";
        globalComm.barrier();
        siloWriter->setDecomposition( 1 );
        siloWriter->writeFile( fname1.str(), 0 );
        globalComm.barrier();
    }
    double t6 = AMP::AMP_MPI::time();

    // Write a seperate output file for each rank
    std::stringstream fname2;
    fname2 << input_file << "_" << globalComm.getSize() << "proc_multiple";
    globalComm.barrier();
    siloWriter->setDecomposition( 2 );
    siloWriter->writeFile( fname2.str(), 0 );
    globalComm.barrier();
    double t7 = AMP::AMP_MPI::time();

    if ( globalComm.getRank() == 0 ) {
        std::cout << "Read in meshes: " << t2 - t1 << std::endl;
        std::cout << "Allocate vectors: " << t3 - t2 << std::endl;
        std::cout << "Register data: " << t4 - t3 << std::endl;
        std::cout << "Initialize vectors: " << t5 - t4 << std::endl;
        if ( globalComm.getSize() <= 20 )
            std::cout << "Write a single file: " << t6 - t5 << std::endl;
        std::cout << "Write multiple files: " << t7 - t6 << std::endl;
        std::cout << "Total time: " << t7 - t1 << std::endl;
    }

    ut->passes( "test ran to completion" );
}


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;
    PROFILE_ENABLE();

#ifdef USE_EXT_SILO
    std::string filename = "input_SiloIO-1";
    if ( argc == 2 )
        filename = argv[1];
    test_Silo( &ut, filename );
#else
    ut.expected_failure( "AMP was not configured with silo" );
#endif

    ut.report();
    PROFILE_SAVE( "test_Silo" );

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
