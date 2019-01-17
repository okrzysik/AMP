#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/utils/shared_ptr.h"
#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <iostream>
#include <string>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );

    // Get the Mesh database and create the mesh parameters
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( params );

    // Construct Variables
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> Variable1(
        new AMP::LinearAlgebra::Variable( "Var1" ) );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> Variable2(
        new AMP::LinearAlgebra::Variable( "Var2" ) );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> Variable3(
        new AMP::LinearAlgebra::Variable( "Var3" ) );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> Variable4(
        new AMP::LinearAlgebra::Variable( "Var4" ) );
    AMP::shared_ptr<AMP::LinearAlgebra::Variable> dummyVar(
        new AMP::LinearAlgebra::Variable( "dummy" ) );

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> subVariable(
        new AMP::LinearAlgebra::MultiVariable( "subVar" ) );
    subVariable->add( Variable1 );
    subVariable->add( Variable2 );

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVariable> fullVariable(
        new AMP::LinearAlgebra::MultiVariable( "fullVariable" ) );
    fullVariable->add( Variable1 );
    fullVariable->add( Variable2 );
    fullVariable->add( Variable3 );
    fullVariable->add( Variable4 );

    // Create the DOF manager
    AMP::Discretization::DOFManager::shared_ptr DOFs =
        AMP::Discretization::simpleDOFManager::create( mesh, AMP::Mesh::GeomType::Vertex, 1, 1 );

    // Create the vectors
    AMP::LinearAlgebra::Vector::shared_ptr multiVector =
        AMP::LinearAlgebra::createVector( DOFs, fullVariable, true );
    AMP::LinearAlgebra::Vector::shared_ptr singleVector =
        AMP::LinearAlgebra::createVector( multiVector->getDOFManager(), dummyVar, false );
    AMP::LinearAlgebra::Vector::shared_ptr subVector =
        multiVector->subsetVectorForVariable( subVariable );
    if ( singleVector->getGlobalSize() == multiVector->getGlobalSize() )
        ut->passes( "single and multivector are the right size" );
    else
        ut->failure( "Sub Vector is the right size" );
    if ( subVector.get() != nullptr ) {
        ut->passes( "Sub Vector is not NULL" );
        if ( multiVector->getGlobalSize() == 2 * subVector->getGlobalSize() )
            ut->passes( "Sub Vector is the right size" );
        else
            ut->failure( "Sub Vector is the right size" );
    } else {
        ut->failure( "Sub Vector is not NULL" );
    }

    // Try to copy data between the single vector and multivector
    singleVector->setRandomValues();
    multiVector->copyVector( singleVector );
    if ( AMP::Utilities::approx_equal( singleVector->L2Norm(), multiVector->L2Norm(), 1e-12 ) )
        ut->passes( "Data copied from single vector to multivector" );
    else
        ut->failure( "Data copied from single vector to multivector" );
    singleVector->zero();
    singleVector->copyVector( multiVector );
    if ( AMP::Utilities::approx_equal( singleVector->L2Norm(), multiVector->L2Norm(), 1e-12 ) )
        ut->passes( "Data copied from multivector to single vector" );
    else
        ut->failure( "Data copied from multivector to single vector" );
}


int testMultiVector( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testMultiVector" );
    for ( auto &exeName : exeNames ) {
        myTest( &ut, exeName );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
