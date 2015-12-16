
#include <utils/AMPManager.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>
#include <utils/Writer.h>

#include <ampmesh/Mesh.h>

#include <vectors/MultiVector.h>
#include <vectors/VectorBuilder.h>
#include <vectors/VectorSelector.h>

#include <discretization/simpleDOF_Manager.h>

#include <operators/map/StridedZAxisMap.h>

#include <functional>


double fooFunctionOfSpace( std::vector<double> const &xyz )
{
    AMP_ASSERT( xyz.size() == 3 );
    return 1.0 + xyz[2];
}

double barFunctionOfSpace( std::vector<double> const &xyz )
{
    AMP_ASSERT( xyz.size() == 3 );
    return 3.0 * std::cos( xyz[2] );
}

void project( AMP::Mesh::MeshIterator const &meshIterator,
              AMP::LinearAlgebra::Vector::shared_ptr vector,
              size_t const dof,
              size_t const dofsPerNode,
              std::function<double( std::vector<double> const & )>
                  functionOfSpace )
{
    AMP_INSIST( dof < dofsPerNode, "WRONG!" );
    AMP::Discretization::DOFManager::const_shared_ptr dofManager = vector->getDOFManager();
    AMP::Mesh::MeshIterator const meshIterator_begin             = meshIterator.begin();
    AMP::Mesh::MeshIterator const meshIterator_end               = meshIterator.end();
    std::vector<size_t> dofIndices;
    std::vector<double> coord;
    for ( AMP::Mesh::MeshIterator iterator = meshIterator_begin; iterator != meshIterator_end;
          ++iterator ) {
        dofManager->getDOFs( iterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == dofsPerNode );
        coord        = iterator->coord();
        double value = functionOfSpace( coord );
        vector->setValueByGlobalID( dofIndices[dof], value );
    } // end for iterator
}


void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string inputFile = "input_" + exeName;
    std::string logFile   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // parse input file
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase( new AMP::InputDatabase( "inputDatabase" ) );
    AMP::InputManager::getManager()->parseInputFile( inputFile, inputDatabase );
    inputDatabase->printClassData( AMP::plog );

    // read the meshe
    AMP::shared_ptr<AMP::Database> meshDatabase = inputDatabase->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( meshDatabase ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh( meshParams );

    AMP::Mesh::Mesh::shared_ptr fooMesh = mesh->Subset( "Foo" );
    AMP::Mesh::Mesh::shared_ptr barMesh = mesh->Subset( "Bar" );

    // build two dof managers
    bool const split     = true;
    int const ghostWidth = 0;
    // TODO: read from input
    int const fooDof         = 3;
    int const barDof         = 0;
    int const fooDofsPerNode = 5;
    int const barDofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr fooDofManager =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::Vertex, ghostWidth, fooDofsPerNode );
    AMP::Discretization::DOFManager::shared_ptr barDofManager =
        AMP::Discretization::simpleDOFManager::create(
            mesh, AMP::Mesh::Vertex, ghostWidth, barDofsPerNode );

    // and two vectors
    AMP::LinearAlgebra::Variable::shared_ptr variable1(
        new AMP::LinearAlgebra::Variable( "fooname" ) ); // TODO: has to match map operator
    AMP::LinearAlgebra::Variable::shared_ptr variable2(
        new AMP::LinearAlgebra::Variable( "barname" ) ); // TODO: has to match map operator
    AMP::LinearAlgebra::Vector::shared_ptr fooFooVector =
        AMP::LinearAlgebra::createVector( fooDofManager, variable1, split );
    AMP::LinearAlgebra::Vector::shared_ptr barFooVector = fooFooVector->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr tmpFooVector = fooFooVector->cloneVector();

    AMP::LinearAlgebra::Vector::shared_ptr barBarVector =
        AMP::LinearAlgebra::createVector( barDofManager, variable2, split );
    AMP::LinearAlgebra::Vector::shared_ptr fooBarVector = barBarVector->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr tmpBarVector = barBarVector->cloneVector();

    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> myVector =
        AMP::LinearAlgebra::MultiVector::create( "MultiSolutionVec", globalComm );
    myVector->addVector( fooFooVector );
    myVector->addVector( barBarVector );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> otherVector =
        AMP::LinearAlgebra::MultiVector::create( "MultiSolutionVec", globalComm );
    otherVector->addVector( tmpFooVector );
    otherVector->addVector( tmpBarVector );

    // fill them
    int const fooBoundaryID = 1;
    int const barBoundaryID = 0;
    AMP::Mesh::MeshIterator fooMeshIterator =
        fooMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, fooBoundaryID );
    project( fooMeshIterator, fooFooVector, fooDof, fooDofsPerNode, fooFunctionOfSpace );
    project( fooMeshIterator, barFooVector, fooDof, fooDofsPerNode, barFunctionOfSpace );

    AMP::Mesh::MeshIterator barMeshIterator =
        barMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, barBoundaryID );
    project( barMeshIterator, barBarVector, barDof, barDofsPerNode, barFunctionOfSpace );
    project( barMeshIterator, fooBarVector, barDof, barDofsPerNode, fooFunctionOfSpace );

    // make map operator
    AMP::shared_ptr<AMP::Database> mapOperatorDatabase =
        inputDatabase->getDatabase( "MapOperator" );
    AMP::shared_ptr<AMP::Operator::Map3to1to3Parameters> mapOperatorParameters(
        new AMP::Operator::Map3to1to3Parameters( mapOperatorDatabase ) );
    mapOperatorParameters->d_Mesh1       = fooMesh;
    mapOperatorParameters->d_Mesh2       = barMesh;
    mapOperatorParameters->d_BoundaryID1 = fooBoundaryID;
    mapOperatorParameters->d_BoundaryID2 = barBoundaryID;
    mapOperatorParameters->d_MapComm     = globalComm;
    AMP::shared_ptr<AMP::Operator::StridedZAxisMap> mapOperator(
        new AMP::Operator::StridedZAxisMap( mapOperatorParameters ) );

    // apply it
    AMP::LinearAlgebra::Vector::shared_ptr dummyVector;

    AMP::pout << "before  "
              << "foo=" << barFooVector->L2Norm() << "  "
              << "bar=" << fooBarVector->L2Norm() << "\n";

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->setDecomposition( 1 );
    siloWriter->registerVector( myVector, mesh, AMP::Mesh::Vertex, "my" );
    siloWriter->registerVector( otherVector, mesh, AMP::Mesh::Vertex, "other" );
    siloWriter->writeFile( "tmp", 0 );
#endif
    mapOperator->setVector( otherVector );

    mapOperator->apply( myVector, dummyVector );

#ifdef USE_EXT_SILO
    siloWriter->writeFile( "tmp", 1 );
#endif

    double const tolerance = 1.0e-14 * otherVector->L2Norm();
    tmpFooVector->subtract( barFooVector, tmpFooVector );
    tmpBarVector->subtract( fooBarVector, tmpBarVector );
    AMP::pout << "after  "
              << "foo=" << tmpFooVector->L2Norm() << "  "
              << "bar=" << tmpBarVector->L2Norm() << "\n";

#ifdef USE_EXT_SILO
    siloWriter->writeFile( "tmp", 2 );
#endif
    double const errNorm = otherVector->L2Norm();
    AMP::pout << "errNorm  " << errNorm << std::endl;

    if ( errNorm < tolerance ) {
        ut->passes( exeName );
    }
    else {
        ut->failure( exeName );
    }
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testStridedZAxisMap" );

    try {
        for ( size_t i = 0; i < exeNames.size(); ++i ) {
            myTest( &ut, exeNames[i] );
        } // end for
    }
    catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    }
    catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
