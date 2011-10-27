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

#include "vectors/Variable.h"
#include "vectors/Vector.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/SiloIO.h"


void test_Silo(AMP::UnitTest *ut) {
    AMP::PIO::logOnlyNodeZero ( "outputMeshManagerTest1" );

#ifdef USE_SILO
    AMP::AMP_MPI commGlobal(AMP_COMM_WORLD);
    commGlobal.barrier();
    int size = commGlobal.getSize();
    double t1 = AMP::AMP_MPI::time();
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( "input_SiloIO" , input_db );
    input_db->printClassData (AMP::plog);

    AMP::Mesh::MeshManagerParameters::shared_ptr  meshmgrParams ( new AMP::Mesh::MeshManagerParameters ( input_db ) );
    AMP::Mesh::MeshManager::shared_ptr  manager ( new AMP::Mesh::MeshManager ( meshmgrParams ) );
    commGlobal.barrier();
    double t2 = AMP::AMP_MPI::time();

    AMP::LinearAlgebra::Variable::shared_ptr  random_var ( new AMP::Mesh::Nodal3VectorVariable ( "random" ) );
    AMP::LinearAlgebra::Variable::shared_ptr  gp_var ( new AMP::Mesh::SingleGaussPointVariable ( "gp_var" ) );
    AMP::LinearAlgebra::Variable::shared_ptr  gp_var2 ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable , 8> ( "gp_var2" ) );
    gp_var->setUnits ( "newton-fathom / acre^2" );
    AMP::LinearAlgebra::Vector::shared_ptr  random = manager->createVector ( random_var );
    AMP::LinearAlgebra::Vector::shared_ptr  gauss_pt = manager->createVector ( gp_var );
    AMP::LinearAlgebra::Vector::shared_ptr  gauss_pt2 = manager->createVector ( gp_var2 );
    AMP::LinearAlgebra::Vector::shared_ptr  displacement = manager->createPositionVector ( "displacement" );
    displacement->getVariable()->setUnits ( "leagues" );
    gauss_pt2->setToScalar ( 100 );

    AMP::LinearAlgebra::Vector::iterator  curd = gauss_pt2->begin();
    AMP::LinearAlgebra::Vector::iterator  end = gauss_pt2->end();
    size_t i = 0;
    while ( curd != end )
    {
        *curd += (double)(i%8);
        i++;
        curd++;
    }

    commGlobal.barrier();
    double t3 = AMP::AMP_MPI::time();
    manager->registerVectorAsData ( displacement );
    manager->registerVectorAsData ( gauss_pt );
    manager->registerVectorAsData ( gauss_pt2 );

    std::stringstream  fname;
    fname << "2pellet_clad_" << size << "proc";
    commGlobal.barrier();
    double t4 = AMP::AMP_MPI::time();
    manager->writeFile<AMP::Mesh::SiloIO> ( fname.str() , 0 );
    commGlobal.barrier();
    double t5 = AMP::AMP_MPI::time();
    for ( int t = 1 ; t <= 2 ; t += 1 )
    {
        gauss_pt->setToScalar ( (double) t );
        random->setRandomValues ();
        displacement->axpby ( .001 , 1.1 , random );
        manager->writeFile<AMP::Mesh::SiloIO> ( fname.str() , 3*t-2 );
        gauss_pt->setToScalar ( 95.9 );
        displacement->setToScalar ( 0.0 );
        manager->writeFile<AMP::Mesh::SiloIO> ( fname.str() , 3*t-1 );
        manager->readFile<AMP::Mesh::SiloIO> ( fname.str() , 3*t-2 );
        manager->writeFile<AMP::Mesh::SiloIO> ( fname.str() , 3*t );
    }
    ut->passes ( "silo file written" );

    commGlobal.barrier();
    double t6 = AMP::AMP_MPI::time();
    int rank = commGlobal.getRank();
    if ( rank == 0 )
    {
        std::cout << "Read in and partition mesh: " << (t2-t1) << std::endl;
        std::cout << "Allocate 2 3-vectors: " << (t3 - t2) << std::endl;
        std::cout << "Write a file: " << (t5-t4) << std::endl;
        std::cout << "10 iterations: " << (t6-t5) << std::endl;
    }
#else
    ut->expected_failure ( "test skipped since no silo configured" );
#endif
}


int main ( int argc , char **argv )
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    test_Silo(&ut);

    ut.report();
    
    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}

