#include <string>
#include <sstream>

#include "utils/Utilities.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/ProfilerApp.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"

#ifdef USE_AMP_VECTORS
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#endif


void test_Silo( AMP::UnitTest *ut, std::string input_file ) {

    AMP::PIO::logOnlyNodeZero ( "output_test_SiloIO" );
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalComm.barrier();
    double t1 = AMP::AMP_MPI::time();

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );
    input_db->printClassData (AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(globalComm);

    // Create the meshes from the input database
    PROFILE_START("Load Mesh");
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(params);
    globalComm.barrier();
    PROFILE_STOP("Load Mesh");
    double t2 = AMP::AMP_MPI::time();

    // Create a surface mesh
    AMP::Mesh::Mesh::shared_ptr submesh = mesh->Subset( mesh->getSurfaceIterator(AMP::Mesh::Face,1) );

#ifdef USE_AMP_VECTORS
    // Create a simple DOFManager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    AMP::Discretization::DOFManager::shared_ptr DOF_scalar = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);
    AMP::Discretization::DOFManager::shared_ptr DOF_vector = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,3,true);
    AMP::Discretization::DOFManager::shared_ptr DOF_gauss  = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Volume,1,8,true);
    AMP::Discretization::DOFManager::shared_ptr DOF_surface = AMP::Discretization::simpleDOFManager::create(submesh,AMP::Mesh::Face,0,1,true);

    // Create the vectors
    AMP::LinearAlgebra::Variable::shared_ptr rank_var( new AMP::LinearAlgebra::Variable("rank") );
    AMP::LinearAlgebra::Vector::shared_ptr rank_vec = AMP::LinearAlgebra::createVector( DOF_scalar, rank_var, true );
    AMP::LinearAlgebra::Variable::shared_ptr position_var( new AMP::LinearAlgebra::Variable("position") );
    AMP::LinearAlgebra::Vector::shared_ptr position = AMP::LinearAlgebra::createVector( DOF_vector, position_var, true );
    AMP::LinearAlgebra::Variable::shared_ptr  gp_var ( new AMP::LinearAlgebra::Variable( "gp_var" ) );
    AMP::LinearAlgebra::Vector::shared_ptr  gauss_pt = AMP::LinearAlgebra::createVector( DOF_gauss, gp_var, true );
    AMP::LinearAlgebra::Variable::shared_ptr  id_var ( new AMP::LinearAlgebra::Variable( "ids" ) );
    AMP::LinearAlgebra::Vector::shared_ptr  id_vec = AMP::LinearAlgebra::createVector( DOF_surface, id_var, true );
    gauss_pt->setToScalar ( 100 );
    globalComm.barrier();
#endif
    double t3 = AMP::AMP_MPI::time();

    // Create a view of a vector
    #ifdef USE_AMP_VECTORS
        AMP::LinearAlgebra::VS_MeshIterator meshSelector( submesh->getIterator(AMP::Mesh::Vertex,1), submesh->getComm() );
        AMP::LinearAlgebra::VS_Stride zSelector(2,3);
        AMP::LinearAlgebra::Vector::shared_ptr  vec_meshSubset = position->select( meshSelector, "mesh subset" );
        AMP_ASSERT(vec_meshSubset.get()!=NULL);
        AMP::LinearAlgebra::Vector::shared_ptr  z_surface = vec_meshSubset->select( zSelector, "z surface" );
        AMP_ASSERT(z_surface.get()!=NULL);
        AMP::Mesh::Mesh::shared_ptr clad = mesh->Subset("clad");
        AMP::LinearAlgebra::Vector::shared_ptr  cladPosition;
        if ( clad.get()!=NULL ) {
            clad->setName("clad");
            AMP::LinearAlgebra::VS_Mesh cladMeshSelector( clad );
            cladPosition = position->select( cladMeshSelector, "cladPosition" );
        }
    #endif

    // Create the silo writer and register the data
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
    int level = 1;  // How much detail do we want to register
    siloWriter->registerMesh( mesh, level );
    siloWriter->registerMesh( submesh, level );
#ifdef USE_AMP_VECTORS
    siloWriter->registerVector( rank_vec, mesh, AMP::Mesh::Vertex, "rank" );
    siloWriter->registerVector( position, mesh, AMP::Mesh::Vertex, "position" );
    siloWriter->registerVector( z_surface, submesh, AMP::Mesh::Vertex, "z_surface" );
    siloWriter->registerVector( gauss_pt, mesh, AMP::Mesh::Volume, "gauss_pnt" );
    siloWriter->registerVector( id_vec, submesh, AMP::Mesh::Face, "surface_ids" );
    // Register a vector over the clad
    if ( clad.get()!=NULL )
        siloWriter->registerVector( cladPosition, clad, AMP::Mesh::Vertex, "clad_position" );
#endif
    globalComm.barrier();
    double t4 = AMP::AMP_MPI::time();


    // Initialize the data
#ifdef USE_AMP_VECTORS
    rank_vec->setToScalar(globalComm.getRank());
    rank_vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    std::vector<size_t> dofs;
    for (AMP::Mesh::MeshIterator it=DOF_vector->getIterator(); it!=it.end(); it++) {
        AMP::Mesh::MeshElementID id = it->globalID();
        DOF_vector->getDOFs( id, dofs );
        std::vector<double> pos = it->coord();
        position->setValuesByGlobalID( dofs.size(), &dofs[0], &pos[0] );
    }
    position->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    id_vec->setToScalar(-1);
    std::vector<int> ids = submesh->getBoundaryIDs();
    for (size_t i=0; i<ids.size(); i++) {
        AMP::Mesh::MeshIterator it = submesh->getBoundaryIDIterator( AMP::Mesh::Face, ids[i], 0 );
        for (size_t j=0; j<it.size(); j++) {
            DOF_surface->getDOFs( it->globalID(), dofs );
            AMP_ASSERT(dofs.size()==1);
            id_vec->setValueByGlobalID( dofs[0], ids[i] );
            ++it;
        }
    }
    globalComm.barrier();
#endif
    double t5 = AMP::AMP_MPI::time();

    // Write a single output file
    if ( globalComm.getSize() <= 20 ) {
        std::stringstream  fname1;
        fname1 << input_file << "_" << globalComm.getSize() << "proc_single";
        globalComm.barrier();
        siloWriter->setDecomposition( 0 );
        siloWriter->writeFile( fname1.str() , 0 );
        globalComm.barrier();
    }
    double t6 = AMP::AMP_MPI::time();

    // Write a seperate output file for each rank
    std::stringstream  fname2;
    fname2 << input_file << "_" << globalComm.getSize() << "proc_multiple";
    globalComm.barrier();
    siloWriter->setDecomposition( 1 );
    siloWriter->writeFile( fname2.str() , 0 );
    globalComm.barrier();
    double t7 = AMP::AMP_MPI::time();

    if ( globalComm.getRank() == 0 ) {
        std::cout << "Read in meshes: " << t2-t1 << std::endl;
        std::cout << "Allocate vectors: " << t3-t2 << std::endl;
        std::cout << "Register data: " << t4-t3 << std::endl;
        std::cout << "Initialize vectors: " << t5-t4 << std::endl;
        if ( globalComm.getSize() <= 20 )
            std::cout << "Write a single file: " << t6-t5 << std::endl;
        std::cout << "Write multiple files: " << t7-t6 << std::endl;
        std::cout << "Total time: " << t7-t1 << std::endl;
    }

    ut->passes("test ran to completion");
}


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    PROFILE_ENABLE();

    #ifdef USE_EXT_SILO
        std::string filename = "input_SiloIO-1";
        if(argc == 2) filename = argv[1];
        test_Silo( &ut, filename );
    #else
        ut.expected_failure("AMP was not configured with silo");
    #endif

    ut.report();
    PROFILE_SAVE("test_Silo");    

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

