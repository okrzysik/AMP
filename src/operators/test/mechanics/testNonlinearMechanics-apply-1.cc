
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <iostream>
#include <string>

#include "discretization/simpleDOF_Manager.h"
#include "vectors/VectorBuilder.h"

#include "materials/Material.h"
#include "operators/LinearOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsNonlinearFEOperator.h"

void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP_INSIST( input_db->keyExists( "testNonlinearMechanicsOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> testNonlinOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db, elementPhysicsModel ) );

    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Variable::shared_ptr var = testNonlinOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( dofMap, var, true );
    AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();

    solVec->setToScalar( 5.0 );

    AMP::pout << "Solution Norm: " << ( solVec->L2Norm() ) << std::endl;

    testNonlinOperator->apply( solVec, resVec );

    double resNorm1 = resVec->L2Norm();

    AMP::pout << "resNorm1 = " << resNorm1 << std::endl;

    testNonlinOperator->apply( solVec, resVec );

    double resNorm2 = resVec->L2Norm();

    AMP::pout << "resNorm2 = " << resNorm2 << std::endl;

    AMP_ASSERT( resNorm1 == resNorm2 );

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNonlinearMechanics-apply-1" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
