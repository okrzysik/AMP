#include <iostream>
#include <string>

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "matrices/MatrixBuilder.h"

#include "operators/OperatorBuilder.h"
#include "operators/OperatorParameters.h"
#include "operators/LinearOperator.h"

void linearTest1( AMP::UnitTest *ut, const std::string &exeName )
{
    // Test create
    const std::string input_file = "input_" + exeName;
    const std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // read the input file into a database
    const auto input_db = AMP::make_shared<AMP::InputDatabase> ( "input_db" );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    // extract the Mesh database and create the mesh parameters
    const auto meshDB = input_db->getDatabase( "Mesh" );
    auto params = AMP::make_shared<AMP::Mesh::MeshParameters>( meshDB );
    params->setComm( globalComm );

    // create the mesh
    const auto meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    // create a linear diffusion operator
    auto linearOperator = AMP::Operator::OperatorBuilder::createOperator( meshAdapter, "LinearDiffusionOp", input_db );
    auto diffOp = AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( linearOperator );

    // concludes creation of a native linear operator
    // ************************************************************************************************
    // extract the internal matrix
    const auto &diffMat = diffOp->getMatrix();

    // extract the left and right vectors
    // COMMENT: these lines will have to be replaced for an external application
    // to provide explicit right and left vectors
    // COMMENT: note that under the hood we are primarily interested in the DOF manager
    // or rather in constructing one. We can't seem to avoid this for an external
    // application
    // COMMENT: We do not give a good description of what left and right vectors are
    // or what information they should provide (elaborate)
    // COMMENT: We do not currently seem to have the ability to construct simple vectors
    // and DOF managers with ghost cells independent of meshes
    const auto leftVector  = diffMat->getLeftVector();
    const auto rightVector = diffMat->getRightVector();

    // COMMENT: this function pointer will need to be set to an actual function
    // COMMENT 2: if we expect users to provide such an interface we should provide one too!!
    std::function<std::vector<size_t>(size_t)> getRow;
    auto newMat = AMP::LinearAlgebra::createMatrix( rightVector, leftVector, "auto", getRow );

    // construct a LinearOperator and set its matrix
    const auto linearOpDB = AMP::make_shared<AMP::InputDatabase> ( "linearOperatorDB" );
    linearOpDB->putInteger("print_info_level", 0);
    auto linearOpParameters = AMP::make_shared<AMP::Operator::OperatorParameters> ( linearOpDB );
    auto linearOp = AMP::make_shared<AMP::Operator::LinearOperator> ( linearOpParameters );
    linearOp->setMatrix( newMat );

    // COMMENT: the next few lines should ideally need to be replaced
    // by a getRowsByGlobalID call that extracts the rows numbered by global ID
    const auto &leftDOFManager =  diffMat->getLeftDOFManager();
    std::vector<unsigned int> uint_cols;
    std::vector<double> values;
    for ( auto row = leftDOFManager->beginDOF(); row < leftDOFManager->endDOF(); ++row ) {
        diffMat->getRowByGlobalID( row, uint_cols, values );
        std::vector<int> cols;
        // COMMENT: for now do an explicit conversion till we fix the interfaces
        cols.assign(uint_cols.begin(), uint_cols.end());  
        // COMMENT: Note the incosistency in the get and set, need to fix!!
        newMat->setValuesByGlobalID( 1, (int) cols.size(), (int *) &row, &cols[0], &values[0]);
    }

    // extract solution, rhs, and residual variables
    const auto uVar = diffOp->getInputVariable();
    const auto v1Var = diffOp->getOutputVariable();
    const auto v2Var = diffOp->getOutputVariable();

    // construct a nodal DOF manager
    const auto nodalScalarDOFManager =
        AMP::Discretization::simpleDOFManager::create( meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    auto  u = AMP::LinearAlgebra::createVector( nodalScalarDOFManager, uVar );
    auto v1 = AMP::LinearAlgebra::createVector( nodalScalarDOFManager, v1Var );
    auto v2 = AMP::LinearAlgebra::createVector( nodalScalarDOFManager, v2Var );

    ut->passes( exeName );

    // Test apply
    bool passed = true;
    for ( int i = 0; i < 10; i++ ) {
        u->setRandomValues();
        v1->setRandomValues();
        v2->setRandomValues();
        diffOp->apply( u, v1 );
        linearOp->apply( u, v2 );
        // COMMENT: simple add, subtract routines would be nice for matrices
        // this test does not really test equivalence, keeping to remind myself
        v2->subtract(v1,v2);
        passed = passed && ( v2->maxNorm() < std::numeric_limits<double>::min());
    } // end for i

    if( passed ) {
        ut->passes( exeName );
    } else {
        ut->failure( "unable to create a copy of a linear operator");
    }
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> files = { "Diffusion-TUI-Thermal-1",     "Diffusion-TUI-Fick-1",
                                       "Diffusion-TUI-Soret-1",       "Diffusion-UO2MSRZC09-Thermal-1",
                                       "Diffusion-UO2MSRZC09-Fick-1", "Diffusion-UO2MSRZC09-Soret-1",
                                       "Diffusion-TUI-TensorFick-1",  "Diffusion-CylindricalFick-1" };

    for ( const auto &file : files )
        linearTest1( &ut, file );

    ut.report();

    const int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
