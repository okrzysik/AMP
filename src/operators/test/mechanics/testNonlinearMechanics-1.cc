
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

void myTest( AMP::UnitTest *ut, const std::string& exeName )
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

    AMP::shared_ptr<AMP::Operator::MechanicsNonlinearFEOperator> testNonlinOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsNonlinearFEOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testNonlinearMechanicsOperator", input_db ) );
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel =
        testNonlinOperator->getMaterialModel();

    AMP_INSIST( input_db->keyExists( "testLinearMechanicsOperator" ), "key missing!" );

    AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperator> testLinOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "testLinearMechanicsOperator", input_db, elementPhysicsModel ) );

    ut->passes( exeName + " : create" );

    AMP::Discretization::DOFManager::shared_ptr dofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );

    AMP::LinearAlgebra::Variable::shared_ptr var = testNonlinOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::createVector( dofMap, var, true );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();

    for ( int j = 0; j < 3; j++ ) {
        solVec->setRandomValues();
        rhsVec->setRandomValues();
        resVec->setRandomValues();
        testNonlinOperator->residual( rhsVec, solVec, resVec );
        resVec->scale( -1.0 );
    } // end for j

    ut->passes( exeName + " : apply" );

    AMP::shared_ptr<AMP::Operator::OperatorParameters> resetParams =
        testNonlinOperator->getParameters( "Jacobian", solVec );

    ut->passes( exeName + " : getJac" );

    testLinOperator->reset( resetParams );

    ut->passes( exeName + " : Linear::reset" );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testNonlinearMechanics-1-normal" );
    exeNames.emplace_back( "testNonlinearMechanics-1-reduced" );
    // exeNames.push_back("testNonlinearMechanics_Creep-1-normal");
    // exeNames.push_back("testNonlinearMechanics_Creep-1-reduced");
    // exeNames.push_back("testNonlinearMechanics_Creep_MatLib-1-normal");
    // exeNames.push_back("testNonlinearMechanics_Creep_MatLib-1-reduced");
    // exeNames.push_back("testNonlinearMechanics_NonlinearStrainHardening-1-normal");
    // exeNames.push_back("testNonlinearMechanics_NonlinearStrainHardening-1-reduced");
    // exeNames.push_back("testVonMises_IsotropicKinematicHardening-1-normal");
    // exeNames.push_back("testNonlinearMechanics_Matpro_Fuel_Creep-1-normal");
    // exeNames.push_back("testFrapconCladCreepMaterialModel-1");
    // exeNames.push_back("testPoroElasticMaterialModel-1");

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
