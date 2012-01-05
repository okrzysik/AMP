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

#include "ampmesh/Mesh.h"
#include "ampmesh/SiloIO.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/NodalVariable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"


void test_Silo( AMP::UnitTest *ut, std::string input_file ) {

    AMP::PIO::logOnlyNodeZero ( "outputMeshManagerTest1" );
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );
    input_db->printClassData (AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(globalComm);

    // Create the meshes from the input database
    boost::shared_ptr<AMP::Mesh::Mesh> mesh = AMP::Mesh::Mesh::buildMesh(params);

    // Create a simple DOFManager
    AMP::Discretization::DOFManagerParameters::shared_ptr DOFparams( new AMP::Discretization::DOFManagerParameters(mesh) );
    AMP::Discretization::DOFManager::shared_ptr DOF_scalar = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,1,true);
    AMP::Discretization::DOFManager::shared_ptr DOF_vector = AMP::Discretization::simpleDOFManager::create(mesh,AMP::Mesh::Vertex,1,3,true);

    // Create the vectors
    AMP::LinearAlgebra::Variable::shared_ptr rank_var( new AMP::Discretization::NodalVariable(1,"rank") );
    AMP::LinearAlgebra::Vector::shared_ptr rank_vec = AMP::LinearAlgebra::createVector( DOF_scalar, rank_var, true );
    AMP::LinearAlgebra::Variable::shared_ptr displacement_var( new AMP::Discretization::NodalVariable(3,"displacement") );
    AMP::LinearAlgebra::Vector::shared_ptr displacement = AMP::LinearAlgebra::createVector( DOF_vector, displacement_var, true );
    //AMP::LinearAlgebra::Variable::shared_ptr  gp_var ( new AMP::Mesh::SingleGaussPointVariable ( "gp_var" ) );
    //AMP::LinearAlgebra::Variable::shared_ptr  gp_var2 ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable , 8> ( "gp_var2" ) );
    //gp_var->setUnits ( "newton-fathom / acre^2" );
    //AMP::LinearAlgebra::Vector::shared_ptr  gauss_pt = manager->createVector ( gp_var );
    //AMP::LinearAlgebra::Vector::shared_ptr  gauss_pt2 = manager->createVector ( gp_var2 );
    //AMP::LinearAlgebra::Vector::shared_ptr  displacement = manager->createPositionVector ( "displacement" );
    //displacement->getVariable()->setUnits ( "leagues" );
    //gauss_pt2->setToScalar ( 100 );

    //AMP::LinearAlgebra::Vector::iterator  curd = gauss_pt2->begin();
    //AMP::LinearAlgebra::Vector::iterator  end = gauss_pt2->end();
    //size_t i = 0;
    //while ( curd != end )
    //{
    //    *curd += (double)(i%8);
    //    i++;
    //    curd++;
    //}

    // Create the silo writer and register the data
    AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
    siloWriter->registerMesh( mesh );
    siloWriter->registerVector( rank_vec, mesh, AMP::Mesh::Vertex, "rank" );
    siloWriter->registerVector( displacement, mesh, AMP::Mesh::Vertex, "displacement" );
    //siloWriter->registerVector( gauss_pt );
    //siloWriter->registerVector( gauss_pt2 );

    // Initialize the data
    rank_vec->setToScalar(globalComm.getRank());
    rank_vec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
    displacement->zero();

    // Write the file
    std::stringstream  fname;
    fname << "2pellet_clad_" << globalComm.getSize() << "proc";
    globalComm.barrier();
    double t4 = AMP::AMP_MPI::time();
    siloWriter->writeFile( fname.str() , 0 );

}


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    #ifdef USE_SILO
        test_Silo( &ut, "input_SiloIO" );
    #else
        ut->expected_failure("AMP was not configured with silo");
    #endif

    ut.report();
    
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

