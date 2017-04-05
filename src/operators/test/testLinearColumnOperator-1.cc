
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include <cstdlib>
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "materials/Material.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearOperator.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/Mesh.h"
#include "discretization/MultiDOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/MultiVariable.h"
#include "vectors/MultiVector.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "applyTests.h"

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testLinearColumnOperator-1" );
    std::string outerInput_file = "input_" + exeName;
    std::string log_file        = "output_" + exeName;
    std::string msgPrefix;

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> outerInput_db( new AMP::InputDatabase( "outerInput_db" ) );
    AMP::InputManager::getManager()->parseInputFile( outerInput_file, outerInput_db );
    outerInput_db->printClassData( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST( outerInput_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> database = outerInput_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    AMP_INSIST( outerInput_db->keyExists( "number_of_tests" ), "key missing!" );
    int numTests = outerInput_db->getInteger( "number_of_tests" );

    for ( int i = 0; i < numTests; i++ ) {
        char key[256];
        sprintf( key, "test_%d", i );

        AMP_INSIST( outerInput_db->keyExists( key ), "key missing!" );
        std::string innerInput_file = outerInput_db->getString( key );
        std::cout << "Running test " << i + 1 << " of " << numTests << ": " << innerInput_file
                  << std::endl;

        AMP::shared_ptr<AMP::InputDatabase> innerInput_db(
            new AMP::InputDatabase( "innerInput_db" ) );
        AMP::InputManager::getManager()->parseInputFile( innerInput_file, innerInput_db );
        innerInput_db->printClassData( AMP::plog );

        AMP_INSIST( innerInput_db->keyExists( "numberOfOperators" ),
                    "key missing!  " + innerInput_file );

        const int numberOfOperators = innerInput_db->getInteger( "numberOfOperators" );

        AMP_INSIST( numberOfOperators >= 1, "more than zero operators need to be specified" );

        AMP_INSIST( innerInput_db->keyExists( "dofsPerNode" ), "key missing!  " + innerInput_file );
        std::vector<int> dofsPerNodeArr = innerInput_db->getIntegerArray( "dofsPerNode" );

        // create a column operator object
        AMP::shared_ptr<AMP::Operator::OperatorParameters> params;
        AMP::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(
            new AMP::Operator::ColumnOperator( params ) );

        std::vector<AMP::LinearAlgebra::Variable::shared_ptr> inputVariables;
        std::vector<AMP::Discretization::DOFManager::shared_ptr> dofMapVec;

        double defTemp = -1.0;
        double defConc = -1.0;
        size_t nVars   = 0;
        for ( int opN = 1; opN <= numberOfOperators; opN++ ) {
            dofMapVec.push_back( AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Vertex, 1, dofsPerNodeArr[opN - 1], true ) );

            char testOpName[256];
            sprintf( testOpName, "testOperator%d", opN );
            AMP_INSIST( innerInput_db->keyExists( testOpName ),
                        "key missing!  " + innerInput_file );

            AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
            AMP::shared_ptr<AMP::Database> testOp_db = innerInput_db->getDatabase( testOpName );
            AMP::shared_ptr<AMP::Operator::Operator> testOperator =
                AMP::Operator::OperatorBuilder::createOperator(
                    meshAdapter, testOpName, innerInput_db, elementPhysicsModel );

            AMP::shared_ptr<AMP::Operator::LinearOperator> myLinOp =
                AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( testOperator );
            AMP_INSIST( myLinOp != nullptr, "Is not a linear operator!" );

            columnOperator->append( testOperator );

            AMP::LinearAlgebra::Variable::shared_ptr opVar = myLinOp->getInputVariable();

            inputVariables.push_back( opVar );

            // this only works as long at least one of the operators is diffusion and
            // its transport model has defaults defined
            AMP::shared_ptr<AMP::Database> model_db;
            if ( testOp_db->keyExists( "VolumeOperator" ) ) {
                AMP::shared_ptr<AMP::Database> volOp_db =
                    innerInput_db->getDatabase( testOp_db->getString( "VolumeOperator" ) );
                if ( ( volOp_db->getName() == "DiffusionNonlinearFEOperator" ) ||
                     ( volOp_db->getName() == "DiffusionLinearFEOperator" ) ) {
                    model_db = innerInput_db->getDatabase( volOp_db->getString( "LocalModel" ) );
                }
            }

            if ( model_db ) {
                defTemp = model_db->getDouble( "Default_Temperature" );
                defConc = model_db->getDouble( "Default_Concentration" );
            }
        }

        msgPrefix = exeName + " : " + innerInput_file;
        ut->passes( msgPrefix + " : create" );

        {
            // Create the vectors
            AMP::LinearAlgebra::Variable::shared_ptr tmp_var(
                new AMP::LinearAlgebra::MultiVariable( "columnInputVariable" ) );
            AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> solVec =
                AMP::LinearAlgebra::MultiVector::create( tmp_var, meshAdapter->getComm() );
            for ( size_t i = 0; i < inputVariables.size(); i++ ) {
                if ( inputVariables[i].get() != nullptr )
                    solVec->addVector(
                        AMP::LinearAlgebra::createVector( dofMapVec[i], inputVariables[i] ) );
            }
            AMP::LinearAlgebra::Vector::shared_ptr rhsVec = solVec->cloneVector();
            AMP::LinearAlgebra::Vector::shared_ptr resVec = solVec->cloneVector();

            for ( size_t i = 0; i < nVars; i++ ) {
                AMP::LinearAlgebra::Variable::shared_ptr opVar = inputVariables[i];
                if ( opVar->getName() == "temperature" ) {
                    AMP::LinearAlgebra::Vector::shared_ptr tVec =
                        solVec->subsetVectorForVariable( inputVariables[i] );
                    tVec->setToScalar( defTemp );
                }
                if ( opVar->getName() == "concentration" ) {
                    AMP::LinearAlgebra::Vector::shared_ptr cVec =
                        solVec->subsetVectorForVariable( inputVariables[i] );
                    cVec->setToScalar( defConc );
                }
            }

            // test apply with single variable vectors
            applyTests( ut, msgPrefix, columnOperator, rhsVec, solVec, resVec );
        }

#if 0
    // test getJacobianParameters
    msgPrefix=exeName + " : " + innerInput_file;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> nullGuess;
    AMP::shared_ptr<AMP::Operator::OperatorParameters> jacobianParameters = testOperator->getParameters("Jacobian", nullGuess);
    if(jacobianParameters.get()!=NULL)
    {
      ut.passes(msgPrefix + "getJacobianParameters (should return NULL for now)");
    }
    else
    {
      ut.numFails++;
    }
#endif

    } // end for i
}

int main( int argc, char *argv[] )
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup( argc, argv, startup_properties );
    AMP::UnitTest ut;

    myTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
