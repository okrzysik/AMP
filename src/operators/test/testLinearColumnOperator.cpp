#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/MultiDOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "applyTests.h"

#include <cstdlib>
#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testLinearColumnOperator-1" );
    std::string outerInput_file = "input_" + exeName;
    std::string log_file        = "output_" + exeName;
    std::string msgPrefix;

    AMP::PIO::logOnlyNodeZero( log_file );

    auto outerInput_db = AMP::Database::parseInputFile( outerInput_file );
    outerInput_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST( outerInput_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto database = outerInput_db->getDatabase( "Mesh" );
    auto params   = std::make_shared<AMP::Mesh::MeshParameters>( database );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    AMP_INSIST( outerInput_db->keyExists( "number_of_tests" ), "key missing!" );
    int numTests = outerInput_db->getScalar<int>( "number_of_tests" );

    for ( int i = 0; i < numTests; i++ ) {
        char key[256];
        sprintf( key, "test_%d", i );

        AMP_INSIST( outerInput_db->keyExists( key ), "key missing!" );
        std::string innerInput_file = outerInput_db->getString( key );
        std::cout << "Running test " << i + 1 << " of " << numTests << ": " << innerInput_file
                  << std::endl;

        auto innerInput_db = AMP::Database::parseInputFile( innerInput_file );
        innerInput_db->print( AMP::plog );

        AMP_INSIST( innerInput_db->keyExists( "numberOfOperators" ),
                    "key missing!  " + innerInput_file );

        const int numberOfOperators = innerInput_db->getScalar<int>( "numberOfOperators" );

        AMP_INSIST( numberOfOperators >= 1, "more than zero operators need to be specified" );

        AMP_INSIST( innerInput_db->keyExists( "dofsPerNode" ), "key missing!  " + innerInput_file );
        std::vector<int> dofsPerNodeArr = innerInput_db->getVector<int>( "dofsPerNode" );

        // create a column operator object
        std::shared_ptr<AMP::Operator::OperatorParameters> colParams;
        auto columnOperator = std::make_shared<AMP::Operator::ColumnOperator>( colParams );

        std::vector<AMP::LinearAlgebra::Variable::shared_ptr> inputVariables;
        std::vector<AMP::Discretization::DOFManager::shared_ptr> dofMapVec;

        double defTemp = -1.0;
        double defConc = -1.0;
        int nVars      = 0;
        for ( int opN = 1; opN <= numberOfOperators; opN++ ) {
            dofMapVec.push_back( AMP::Discretization::simpleDOFManager::create(
                meshAdapter, AMP::Mesh::GeomType::Vertex, 1, dofsPerNodeArr[opN - 1], true ) );

            char testOpName[256];
            sprintf( testOpName, "testOperator%d", opN );
            AMP_INSIST( innerInput_db->keyExists( testOpName ),
                        "key missing!  " + innerInput_file );

            std::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
            auto testOp_db    = innerInput_db->getDatabase( testOpName );
            auto testOperator = AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, testOpName, innerInput_db, elementPhysicsModel );

            auto myLinOp = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( testOperator );
            AMP_INSIST( myLinOp != nullptr, "Is not a linear operator!" );

            columnOperator->append( testOperator );

            AMP::LinearAlgebra::Variable::shared_ptr opVar = myLinOp->getInputVariable();

            inputVariables.push_back( opVar );

            // this only works as long at least one of the operators is diffusion and
            // its transport model has defaults defined
            std::shared_ptr<AMP::Database> model_db;
            if ( testOp_db->keyExists( "VolumeOperator" ) ) {
                auto volOp_db =
                    innerInput_db->getDatabase( testOp_db->getString( "VolumeOperator" ) );
                if ( ( volOp_db->getName() == "DiffusionNonlinearFEOperator" ) ||
                     ( volOp_db->getName() == "DiffusionLinearFEOperator" ) ) {
                    model_db = innerInput_db->getDatabase( volOp_db->getString( "LocalModel" ) );
                }
            }

            if ( model_db ) {
                defTemp = model_db->getScalar<double>( "Default_Temperature" );
                defConc = model_db->getScalar<double>( "Default_Concentration" );
            }
        }

        msgPrefix = exeName + " : " + innerInput_file;
        ut->passes( msgPrefix + " : create" );

        {
            // Create the vectors
            auto tmp_var =
                std::make_shared<AMP::LinearAlgebra::MultiVariable>( "columnInputVariable" );
            std::shared_ptr<AMP::LinearAlgebra::MultiVector> solVec =
                AMP::LinearAlgebra::MultiVector::create( tmp_var, meshAdapter->getComm() );
            for ( size_t iv = 0; iv < inputVariables.size(); iv++ ) {
                if ( inputVariables[iv].get() != nullptr )
                    solVec->addVector(
                        AMP::LinearAlgebra::createVector( dofMapVec[iv], inputVariables[iv] ) );
            }
            auto rhsVec = solVec->cloneVector();
            auto resVec = solVec->cloneVector();

            for ( int iv = 0; iv < nVars; iv++ ) {
                auto opVar = inputVariables[iv];
                if ( opVar->getName() == "temperature" ) {
                    auto tVec = solVec->subsetVectorForVariable( inputVariables[iv] );
                    tVec->setToScalar( defTemp, tVec );
                }
                if ( opVar->getName() == "concentration" ) {
                    auto cVec = solVec->subsetVectorForVariable( inputVariables[iv] );
                    cVec->setToScalar( defConc, cVec );
                }
            }

            // test apply with single variable vectors
            applyTests( ut, msgPrefix, columnOperator, rhsVec, solVec, resVec );
        }

#if 0
    // test getJacobianParameters
    msgPrefix=exeName + " : " + innerInput_file;
    std::shared_ptr<AMP::LinearAlgebra::Vector> nullGuess;
    std::shared_ptr<AMP::Operator::OperatorParameters> jacobianParameters = testOperator->getParameters("Jacobian", nullGuess);
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

int testLinearColumnOperator( int argc, char *argv[] )
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
