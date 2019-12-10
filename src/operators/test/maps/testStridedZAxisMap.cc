#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/map/StridedZAxisMap.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"

#include <functional>


double fooFunctionOfSpace( AMP::Mesh::Point const &xyz )
{
    AMP_ASSERT( xyz.size() == 3 );
    return 1.0 + xyz[2];
}

double barFunctionOfSpace( AMP::Mesh::Point const &xyz )
{
    AMP_ASSERT( xyz.size() == 3 );
    return 3.0 * std::cos( xyz[2] );
}

static void project( AMP::Mesh::MeshIterator const &meshIterator,
                     AMP::LinearAlgebra::Vector::shared_ptr vector,
                     size_t const dof,
                     size_t const dofsPerNode,
                     std::function<double( AMP::Mesh::Point const & )> functionOfSpace )
{
    AMP_INSIST( dof < dofsPerNode, "WRONG!" );
    auto dofManager         = vector->getDOFManager();
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    std::vector<size_t> dofIndices;
    for ( auto iterator = meshIterator_begin; iterator != meshIterator_end; ++iterator ) {
        dofManager->getDOFs( iterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == dofsPerNode );
        auto coord   = iterator->coord();
        double value = functionOfSpace( coord );
        vector->setValueByGlobalID( dofIndices[dof], value );
    } // end for iterator
}


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string inputFile = "input_" + exeName;
    std::string logFile   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( logFile );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // parse input file
    auto inputDatabase = AMP::Database::parseInputFile( inputFile );
    inputDatabase->print( AMP::plog );

    // read the meshe
    auto meshDatabase = inputDatabase->getDatabase( "Mesh" );
    auto meshParams   = AMP::make_shared<AMP::Mesh::MeshParameters>( meshDatabase );
    meshParams->setComm( globalComm );
    auto mesh    = AMP::Mesh::Mesh::buildMesh( meshParams );
    auto fooMesh = mesh->Subset( "Foo" );
    auto barMesh = mesh->Subset( "Bar" );

    // build two dof managers
    bool const split     = true;
    int const ghostWidth = 0;
    // TODO: read from input
    int const fooDof         = 3;
    int const barDof         = 0;
    int const fooDofsPerNode = 5;
    int const barDofsPerNode = 1;
    auto fooDofManager       = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, ghostWidth, fooDofsPerNode );
    auto barDofManager = AMP::Discretization::simpleDOFManager::create(
        mesh, AMP::Mesh::GeomType::Vertex, ghostWidth, barDofsPerNode );

    // and two vectors
    auto variable1    = AMP::make_shared<AMP::LinearAlgebra::Variable>( "fooname" );
    auto variable2    = AMP::make_shared<AMP::LinearAlgebra::Variable>( "barname" );
    auto fooFooVector = AMP::LinearAlgebra::createVector( fooDofManager, variable1, split );
    auto barFooVector = fooFooVector->cloneVector();
    auto tmpFooVector = fooFooVector->cloneVector();

    auto barBarVector = AMP::LinearAlgebra::createVector( barDofManager, variable2, split );
    auto fooBarVector = barBarVector->cloneVector();
    auto tmpBarVector = barBarVector->cloneVector();

    auto myVector = AMP::LinearAlgebra::MultiVector::create( "MultiSolutionVec", globalComm );
    myVector->addVector( fooFooVector );
    myVector->addVector( barBarVector );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> otherVector =
        AMP::LinearAlgebra::MultiVector::create( "MultiSolutionVec", globalComm );
    otherVector->addVector( tmpFooVector );
    otherVector->addVector( tmpBarVector );

    // fill them
    int const fooBoundaryID = 1;
    int const barBoundaryID = 0;
    auto fooMeshIterator =
        fooMesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, fooBoundaryID );
    project( fooMeshIterator, fooFooVector, fooDof, fooDofsPerNode, fooFunctionOfSpace );
    project( fooMeshIterator, barFooVector, fooDof, fooDofsPerNode, barFunctionOfSpace );

    auto barMeshIterator =
        barMesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, barBoundaryID );
    project( barMeshIterator, barBarVector, barDof, barDofsPerNode, barFunctionOfSpace );
    project( barMeshIterator, fooBarVector, barDof, barDofsPerNode, fooFunctionOfSpace );

    // make map operator
    auto mapOperatorDatabase = inputDatabase->getDatabase( "MapOperator" );
    auto mapOperatorParameters =
        AMP::make_shared<AMP::Operator::Map3to1to3Parameters>( mapOperatorDatabase );
    mapOperatorParameters->d_Mesh1       = fooMesh;
    mapOperatorParameters->d_Mesh2       = barMesh;
    mapOperatorParameters->d_BoundaryID1 = fooBoundaryID;
    mapOperatorParameters->d_BoundaryID2 = barBoundaryID;
    mapOperatorParameters->d_MapComm     = globalComm;
    auto mapOperator = AMP::make_shared<AMP::Operator::StridedZAxisMap>( mapOperatorParameters );

    // apply it
    AMP::LinearAlgebra::Vector::shared_ptr dummyVector;

    AMP::pout << "before  "
              << "foo=" << barFooVector->L2Norm() << "  "
              << "bar=" << fooBarVector->L2Norm() << "\n";

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->setDecomposition( 1 );
    siloWriter->registerVector( myVector, mesh, AMP::Mesh::GeomType::Vertex, "my" );
    siloWriter->registerVector( otherVector, mesh, AMP::Mesh::GeomType::Vertex, "other" );
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
    } else {
        ut->failure( exeName );
    }
}


int testStridedZAxisMap( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testStridedZAxisMap" );

    for ( auto &exeName : exeNames )
        myTest( &ut, exeName );

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
