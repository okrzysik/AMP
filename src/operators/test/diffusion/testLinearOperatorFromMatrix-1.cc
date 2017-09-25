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
#include "vectors/SimpleVector.h"
#include "vectors/VectorBuilder.h"

#include "matrices/MatrixBuilder.h"

#include "operators/OperatorBuilder.h"
#include "operators/OperatorParameters.h"
#include "operators/LinearOperator.h"

void userLinearOperatorTest( AMP::UnitTest * const ut, const std::string &exeName )
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
    
    // extract the internal matrix
    const auto &userMat = diffOp->getMatrix();

    AMP_INSIST( userMat->numGlobalColumns()==userMat->numGlobalRows(), "matrix is not square" );

    // extract the right vector
    const auto userVector  = userMat->getRightVector();

    // concludes creation of a native linear operator
    // ************************************************************************************************
    
    // extract information about the local size and mpi comm
    const auto localSize    = userVector->getLocalSize();
    const auto ampComm      = userVector->getComm();

    // construct a dof manager
    const auto dofManager = AMP::make_shared<AMP::Discretization::DOFManager> ( localSize, ampComm );  
    const auto copyVariable = AMP::make_shared<AMP::LinearAlgebra::Variable> ( "copyVariable" ); 
    
    // create a vector based on the dofs and variable
    auto ampVector = AMP::LinearAlgebra::createVector( dofManager, copyVariable );
    AMP_INSIST(ampVector!=nullptr, "ampVector is null" );

    // copy values from one vector to another    
    std::copy( userVector->begin(), userVector->end(), ampVector->begin());
    ampVector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    // concludes demonstrating how to initialize an AMP vector from a user vector
    // ************************************************************************************************
 
    // create a lambda that returns non zero column ids given a global row id
    auto getColumnIDS = [ userMat ]( size_t row ) { return userMat->getColumnIDs(row); };

    // create a matrix based on the dimensions of the copied vector
    auto ampMat = AMP::LinearAlgebra::createMatrix( ampVector, ampVector, "auto",  getColumnIDS );

    // construct a LinearOperator and set its matrix
    const auto linearOpDB = AMP::make_shared<AMP::InputDatabase> ( "linearOperatorDB" );
    linearOpDB->putInteger("print_info_level", 0);
    auto linearOpParameters = AMP::make_shared<AMP::Operator::OperatorParameters> ( linearOpDB );
    auto linearOp = AMP::make_shared<AMP::Operator::LinearOperator> ( linearOpParameters );
    linearOp->setMatrix( ampMat );
    linearOp->setVariables( copyVariable, copyVariable );

    // copy the user matrix into the amp matrix
    std::vector<double> coefficients;
    std::vector<size_t> cols;
    const size_t numRows = 1;
    for ( auto row = userMat->beginRow(); row < userMat->endRow(); ++row ) {
        userMat->getRowByGlobalID( row, cols, coefficients );        
        ampMat->setValuesByGlobalID( numRows, cols.size(), &row, cols.data(), coefficients.data());
    }
    // concludes demonstrating how to initialize an AMP linear operator from a user matrix
    // ************************************************************************************************
    
    auto u = ampVector->cloneVector();
    auto v = ampVector->cloneVector();

    u->setRandomValues();
    v->setRandomValues();

    // form the difference of the matrices
    // COMMENT: simple add, subtract routines would be nice for matrices
    ampMat->axpy(-1.0, userMat);
    
    linearOp->apply( u, v );

    auto passed = ( v->maxNorm() <= std::numeric_limits<double>::min());

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

    std::vector<std::string> files = { "Diffusion-TUI-Thermal-1",     "Diffusion-UO2MSRZC09-Thermal-1" };

    for ( const auto &file : files )
        userLinearOperatorTest( &ut, file );

    ut.report();

    const int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
