
#include <utils/AMPManager.h>
#include <utils/UnitTest.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/Writer.h>
#include <utils/PIO.h>

#include <ampmesh/Mesh.h>

#include <vectors/VectorBuilder.h>
#include <vectors/VectorSelector.h>
#include <vectors/MultiVector.h>

#include <discretization/simpleDOF_Manager.h>

#include <operators/StridedZAxisMap.h>

#include <boost/function.hpp>

double fooFunctionOfSpace(std::vector<double> const & xyz) 
{ 
    AMP_ASSERT( xyz.size() == 3 );
    return 1.0 + xyz[2]; 
}

double barFunctionOfSpace(std::vector<double> const & xyz) 
{ 
    AMP_ASSERT( xyz.size() == 3 );
    return 3.0 * std::cos(xyz[2]); 
}

void project
    ( AMP::Mesh::MeshIterator const &                      meshIterator
    , AMP::LinearAlgebra::Vector::shared_ptr               vector
    , size_t const                                         dof
    , size_t const                                         dofsPerNode
    , boost::function<double(std::vector<double> const &)> functionOfSpace
    )
{
    AMP_INSIST( dof < dofsPerNode, "WRONG!" );
    AMP::Discretization::DOFManager::const_shared_ptr dofManager = vector->getDOFManager();
    AMP::Mesh::MeshIterator const meshIterator_begin = meshIterator.begin();
    AMP::Mesh::MeshIterator const meshIterator_end   = meshIterator.end();
    std::vector<size_t> dofIndices;
    std::vector<double> coord;
    double value;
    for (AMP::Mesh::MeshIterator iterator = meshIterator_begin; 
        iterator != meshIterator_end;
        ++iterator) {
        dofManager->getDOFs(iterator->globalID(), dofIndices);
        AMP_ASSERT( dofIndices.size() == dofsPerNode );
        coord = iterator->coord();
        value = functionOfSpace(coord);
        vector->setValueByGlobalID(dofIndices[dof], value);
    } // end for iterator
}
    

void myTest(AMP::UnitTest *ut, std::string exeName) 
{
    std::string inputFile = "input_" + exeName;
    std::string logFile = "output_" + exeName; 
  
    AMP::PIO::logOnlyNodeZero(logFile);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  
    // parse input file
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase(new AMP::InputDatabase("inputDatabase"));
    AMP::InputManager::getManager()->parseInputFile(inputFile, inputDatabase);
    inputDatabase->printClassData(AMP::plog);
  
    // read the meshe
    AMP::shared_ptr<AMP::Database> meshDatabase = inputDatabase->getDatabase("Mesh");
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(globalComm);
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);
  
    AMP::Mesh::Mesh::shared_ptr fooMesh = mesh->Subset("Foo");
    AMP::Mesh::Mesh::shared_ptr barMesh = mesh->Subset("Bar");

    // build two dof managers
    bool const split = true;
    int const ghostWidth = 0;
    // TODO: read from input
    int const fooDof = 3;
    int const barDof = 0;
    int const fooDofsPerNode = 5;
    int const barDofsPerNode = 1;
    AMP::Discretization::DOFManager::shared_ptr fooDofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, fooDofsPerNode);
    AMP::Discretization::DOFManager::shared_ptr barDofManager = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, barDofsPerNode);

    // and two vectors
    AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("noname")); // TODO: has to match map operator
    AMP::LinearAlgebra::Vector::shared_ptr fooVector = AMP::LinearAlgebra::createVector(fooDofManager, variable, split);
    AMP::LinearAlgebra::Vector::shared_ptr barVector = AMP::LinearAlgebra::createVector(barDofManager, variable, split);
    
    // fill them
    int const fooBoundaryID = 1;
    int const barBoundaryID = 0;
    AMP::Mesh::MeshIterator fooMeshIterator = fooMesh->getBoundaryIDIterator(AMP::Mesh::Vertex, fooBoundaryID);
    project(fooMeshIterator, fooVector, fooDof, fooDofsPerNode, fooFunctionOfSpace);

    AMP::Mesh::MeshIterator barMeshIterator = barMesh->getBoundaryIDIterator(AMP::Mesh::Vertex, barBoundaryID);
    project(barMeshIterator, barVector, barDof, barDofsPerNode, barFunctionOfSpace);

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
    siloWriter->setDecomposition(1);
    siloWriter->registerVector(fooVector, fooMesh, AMP::Mesh::Vertex, "foo");
    siloWriter->registerVector(barVector, barMesh, AMP::Mesh::Vertex, "bar");
    siloWriter->writeFile("tmp", 0);
#endif

    AMP_INSIST( barDofsPerNode == 1, "We'll see about that..." );
    AMP::LinearAlgebra::Vector::shared_ptr errorVector = barVector->cloneVector();
    errorVector->subtract(barVector, fooVector->constSelect(AMP::LinearAlgebra::VS_Stride(fooDof, fooDofsPerNode), variable->getName()));
    AMP::pout<<"before="<<errorVector->L2Norm()<<"\n";

    // make map operator
    AMP::shared_ptr<AMP::Database> mapOperatorDatabase = inputDatabase->getDatabase("MapOperator");
    AMP::shared_ptr<AMP::Operator::Map3to1to3Parameters> mapOperatorParameters(new AMP::Operator::Map3to1to3Parameters(mapOperatorDatabase));
    mapOperatorParameters->d_Mesh1= fooMesh;
    mapOperatorParameters->d_Mesh2= barMesh;
    mapOperatorParameters->d_BoundaryID1= fooBoundaryID;
    mapOperatorParameters->d_BoundaryID2= barBoundaryID;
    mapOperatorParameters->d_MapComm = globalComm;
    AMP::shared_ptr<AMP::Operator::StridedZAxisMap> mapOperator(new AMP::Operator::StridedZAxisMap(mapOperatorParameters));

    // apply it
    AMP::LinearAlgebra::Vector::shared_ptr dummyVector;
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> tmpVector;
    tmpVector.push_back(fooVector);
    tmpVector.push_back(barVector);
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> fooBarVector = AMP::LinearAlgebra::MultiVector::create(variable, globalComm, tmpVector);
//    AMP::LinearAlgebra::MultiVector::shared_ptr fooBarVector = AMP::LinearAlgebra::MultiVector::create(variable, globalComm);
//    fooBarVector->encapsulate(fooVector);
//    fooBarVector->encapsulate(barVector);
    AMP::LinearAlgebra::Vector::shared_ptr otherVector = fooBarVector->cloneVector();
    mapOperator->setVector(otherVector);
    
    mapOperator->apply(dummyVector, fooBarVector, dummyVector);

    errorVector->subtract(barVector, fooVector->constSelect(AMP::LinearAlgebra::VS_Stride(fooDof, fooDofsPerNode), variable->getName()));
    AMP::pout<<"after="<<errorVector->L2Norm()<<"\n";

    double const tolerance = 1.0e-14 * barVector->L2Norm();
    if (errorVector->L2Norm() < tolerance) {
        ut->passes(exeName);
    } else {
        ut->failure(exeName);
    }
}


int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
  AMP::UnitTest ut;

  std::vector<std::string> exeNames; 
  exeNames.push_back("testStridedZAxisMap");

  try {
    for (size_t i = 0; i < exeNames.size(); ++i) { 
      myTest(&ut, exeNames[i]); 
    } // end for
  } catch (std::exception &err) {
    std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();
  int num_failed = ut.NumFailGlobal();

  AMP::AMPManager::shutdown();
  return num_failed;
}  



