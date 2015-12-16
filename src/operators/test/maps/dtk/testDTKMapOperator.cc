
#include <ampmesh/Mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <utils/AMPManager.h>
#include <utils/AMP_MPI.h>
#include <utils/Database.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>
#include <utils/Utilities.h>
#include <utils/Writer.h>
#include <utils/shared_ptr.h>
#include <vectors/VectorBuilder.h>

#include <operators/map/dtk/DTKMapOperator.h>

#include <cstdlib>
#include <iostream>
#include <string>

double testFunction1( const std::vector<double> &coords )
{
    return coords[0] + coords[1] + coords[2] + 2.0;
}

void myTest( AMP::UnitTest *ut )
{
    std::string exeName( "testDTKMapOperator" );
    std::string log_file = "output_" + exeName;
    std::string msgPrefix;
    AMP::PIO::logOnlyNodeZero( log_file );

    // load the source mesh
    AMP::pout << "Loading the source mesh" << std::endl;
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    std::string input_file = "input_" + exeName;
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::shared_ptr<AMP::Database> sourceMeshDatabase = input_db->getDatabase( "SourceMesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> sourceMeshParams(
        new AMP::Mesh::MeshParameters( sourceMeshDatabase ) );
    sourceMeshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr sourceMesh    = AMP::Mesh::Mesh::buildMesh( sourceMeshParams );
    std::size_t const numVerticesOnSourceMesh = sourceMesh->numGlobalElements( AMP::Mesh::Vertex );
    std::size_t const numElementsOnSourceMesh = sourceMesh->numGlobalElements( AMP::Mesh::Volume );
    AMP::pout << "source mesh contains " << numVerticesOnSourceMesh << " vertices\n";
    AMP::pout << "source mesh contains " << numElementsOnSourceMesh << " elements\n";

    // build source vector
    AMP::pout << "Building the source vector" << std::endl;
    bool const split      = true;
    int const ghostWidth  = 1;
    int const dofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr sourceDofManager =
        AMP::Discretization::simpleDOFManager::create(
            sourceMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode );
    AMP::LinearAlgebra::Variable::shared_ptr variable(
        new AMP::LinearAlgebra::Variable( "dummy" ) );
    AMP::LinearAlgebra::Vector::shared_ptr sourceVector =
        AMP::LinearAlgebra::createVector( sourceDofManager, variable, split );
    // and fill it
    std::vector<std::size_t> dofIndices;
    double value;
    AMP::pout << "Filling source vector" << std::endl;
    AMP::Mesh::MeshIterator sourceMeshIterator = sourceMesh->getIterator( AMP::Mesh::Vertex );
    for ( sourceMeshIterator = sourceMeshIterator.begin();
          sourceMeshIterator != sourceMeshIterator.end();
          ++sourceMeshIterator ) {
        sourceDofManager->getDOFs( sourceMeshIterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        value = testFunction1( sourceMeshIterator->coord() );
        sourceVector->setLocalValueByGlobalID( dofIndices[0], value );
    }
#ifdef USE_EXT_SILO
    {
        AMP::Utilities::Writer::shared_ptr siloWriter =
            AMP::Utilities::Writer::buildWriter( "silo" );
        siloWriter->setDecomposition( 1 );
        siloWriter->registerVector( sourceVector, sourceMesh, AMP::Mesh::Vertex, "vector" );
        siloWriter->writeFile( "source", 0 );
    }
#endif

    // load the target mesh
    AMP::pout << "Loading the target mesh" << std::endl;
    AMP::shared_ptr<AMP::Database> targetMeshDatabase = input_db->getDatabase( "TargetMesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> targetMeshParams(
        new AMP::Mesh::MeshParameters( targetMeshDatabase ) );
    targetMeshParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    AMP::Mesh::Mesh::shared_ptr targetMesh    = AMP::Mesh::Mesh::buildMesh( targetMeshParams );
    std::size_t const numVerticesOnTargetMesh = targetMesh->numGlobalElements( AMP::Mesh::Vertex );
    std::size_t const numElementsOnTargetMesh = targetMesh->numGlobalElements( AMP::Mesh::Volume );
    AMP::pout << "target mesh contains " << numVerticesOnTargetMesh << " vertices\n";
    AMP::pout << "target mesh contains " << numElementsOnTargetMesh << " elements\n";

    AMP::pout << "Building the target vector" << std::endl;
    AMP::Discretization::DOFManager::shared_ptr targetDofManager =
        AMP::Discretization::simpleDOFManager::create(
            targetMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode );
    AMP::LinearAlgebra::Vector::shared_ptr targetVector =
        AMP::LinearAlgebra::createVector( targetDofManager, variable, split );

    // create dtk map operator.
    AMP::shared_ptr<AMP::Database> null_db;
    AMP::shared_ptr<AMP::Operator::DTKMapOperatorParameters> dtk_op_params(
        new AMP::Operator::DTKMapOperatorParameters( null_db ) );
    dtk_op_params->d_domain_mesh = sourceMesh;
    dtk_op_params->d_range_mesh  = targetMesh;
    dtk_op_params->d_domain_dofs = sourceDofManager;
    dtk_op_params->d_range_dofs  = targetDofManager;
    AMP::shared_ptr<AMP::Operator::Operator> dtk_operator(
        new AMP::Operator::DTKMapOperator( dtk_op_params ) );

    // apply the map.
    AMP::pout << "Apply dtk operator" << std::endl;
    AMP::LinearAlgebra::Vector::shared_ptr null_vector;
    dtk_operator->apply( null_vector, sourceVector, targetVector );

    // checking the answer
    AMP::pout << "Check answer" << std::endl;
    AMP::pout << "source vector l2 norm = " << sourceVector->L2Norm() << std::endl;
    AMP::pout << "target vector l2 norm = " << targetVector->L2Norm() << std::endl;
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "silo" );
    siloWriter->setDecomposition( 1 );
    siloWriter->registerVector( targetVector, targetMesh, AMP::Mesh::Vertex, "vector" );
    siloWriter->writeFile( "target", 0 );
#endif
    double const atol                          = 1.0e-14;
    double const rtol                          = 1.0e-14;
    double const tol                           = atol + rtol * targetVector->L2Norm();
    AMP::Mesh::MeshIterator targetMeshIterator = targetMesh->getIterator( AMP::Mesh::Vertex );
    for ( targetMeshIterator = targetMeshIterator.begin();
          targetMeshIterator != targetMeshIterator.end();
          ++targetMeshIterator ) {
        targetDofManager->getDOFs( targetMeshIterator->globalID(), dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        value = testFunction1( targetMeshIterator->coord() );
        targetVector->addLocalValueByGlobalID( dofIndices[0], -value );
    }
    AMP::pout << "error l2 norm = " << targetVector->L2Norm() << std::endl;
#ifdef USE_EXT_SILO
    siloWriter->writeFile( "target", 1 );
#endif

    AMP_ASSERT( targetVector->L2Norm() < tol );

    ut->passes( exeName );
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
