#include <utils/AMPManager.h>
#include <utils/UnitTest.h>
#include <utils/InputDatabase.h>
#include <utils/Writer.h>
#include <utils/PIO.h>
#include <ampmesh/mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <vectors/VectorBuilder.h>
#include <vectors/MultiVariable.h>
#include <vectors/VectorSelector.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <DTK_MeshFreeInterpolator.hpp>
#include <DTK_MeshFreeInterpolatorFactory.hpp>

void getVerticesCoordinates(
    AMP::Mesh::MeshIterator                     meshIterator
  , AMP::Discretization::DOFManager::shared_ptr dofManager
  , AMP::LinearAlgebra::Vector::shared_ptr      xyzVector
  )
{
  size_t const numLocalVertices = meshIterator.size();
  std::vector<size_t> dofIndices;
  std::vector<double> coordinates;
  AMP::Mesh::MeshIterator meshIterator_begin = meshIterator.begin();
  AMP::Mesh::MeshIterator meshIterator_end = meshIterator.end();
  AMP_ASSERT(meshIterator == meshIterator_begin);
  for (size_t vertex = 0; vertex < numLocalVertices; ++vertex, ++meshIterator) {
    coordinates = meshIterator->coord();
    dofManager->getDOFs(meshIterator->globalID(), dofIndices);
    AMP_ASSERT(dofIndices.size() == 3);
    xyzVector->setLocalValuesByGlobalID(3, &(dofIndices[0]), &(coordinates[0]));
  }
  AMP_ASSERT(meshIterator == meshIterator_end);
}

void rearrangeData(
    AMP::LinearAlgebra::Vector::const_shared_ptr stridedVector
  , size_t const                                 stride
  , std::vector<double> &                        buffer
  )
{
  size_t const numLocalVertices = stridedVector->getLocalSize() / stride;
  buffer.resize(stride*numLocalVertices);
  AMP_ASSERT(stridedVector->getLocalSize() == buffer.size());
  AMP::LinearAlgebra::Vector::const_shared_ptr tmpVector;
  for (size_t i = 0; i < stride; ++i) {
    tmpVector = stridedVector->constSelect(AMP::LinearAlgebra::VS_Stride(i, stride), "tmp");
    AMP_ASSERT(tmpVector->getLocalSize() == numLocalVertices);
    tmpVector->copyOutRawData(&(buffer[0])+i*numLocalVertices);
  }
}

void rearrangeData(
    std::vector<double> const &            buffer
  , size_t const                           stride
  , AMP::LinearAlgebra::Vector::shared_ptr stridedVector
  )
{
  size_t const numLocalVertices = buffer.size() / stride;
  AMP_ASSERT(stridedVector->getLocalSize() == buffer.size());
  AMP::LinearAlgebra::Vector::shared_ptr tmpVector;
  for (size_t i = 0; i < stride; ++i) {
    tmpVector = stridedVector->select(AMP::LinearAlgebra::VS_Stride(i, stride), "tmp");
    AMP_ASSERT(tmpVector->getLocalSize() == numLocalVertices);
    tmpVector->putRawData(&(buffer[0])+i*numLocalVertices);
  }
}

void getData(
    AMP::Mesh::MeshIterator                     meshIterator
  , AMP::Discretization::DOFManager::shared_ptr dofManager
  , std::vector<double> &                       buffer
  )
{
  bool const split = true;
  AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("xyz"));
  AMP::LinearAlgebra::Vector::shared_ptr xyzVector = AMP::LinearAlgebra::createVector(dofManager, variable, split);
  getVerticesCoordinates(meshIterator, dofManager, xyzVector);
  buffer.resize(xyzVector->getLocalSize());
  xyzVector->copyOutRawData(&(buffer[0]));
}

void setupDTKInterpolator(
    AMP::Mesh::Mesh::shared_ptr                         sourceMesh
  , AMP::Mesh::Mesh::shared_ptr                         targetMesh
  , Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator
  , double const                                        radius
  )
{
  // TODO: do we want to pass the dof managers as argument?
  int const ghostWidth = 0;
  int const dofsPerNode = 3;
  AMP::Discretization::DOFManager::shared_ptr sourceDofManager = AMP::Discretization::simpleDOFManager::create(sourceMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
  AMP::Discretization::DOFManager::shared_ptr targetDofManager = AMP::Discretization::simpleDOFManager::create(targetMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);

  std::vector<double> sourceBuffer;
  std::vector<double> targetBuffer;
  getData(sourceMesh->getIterator(AMP::Mesh::Vertex, ghostWidth), sourceDofManager, sourceBuffer);
  getData(targetMesh->getIterator(AMP::Mesh::Vertex, ghostWidth), targetDofManager, targetBuffer);
  Teuchos::ArrayView<double> sourcePoints(&(sourceBuffer[0]), sourceBuffer.size());
  Teuchos::ArrayView<double> targetPoints(&(targetBuffer[0]), targetBuffer.size());
  interpolator->setProblem(sourcePoints, targetPoints, radius);
}

void applyDTKInterpolator(
    Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator
  , AMP::LinearAlgebra::Vector::const_shared_ptr        sourceVector
  , AMP::LinearAlgebra::Vector::shared_ptr              targetVector
  , size_t const                                        data_dim
  )
{
  std::vector<double> sourceBuffer(sourceVector->getLocalSize());
  std::vector<double> targetBuffer(targetVector->getLocalSize());
  Teuchos::ArrayView<double> sourceData(&(sourceBuffer[0]), sourceBuffer.size());
  Teuchos::ArrayView<double> targetData(&(targetBuffer[0]), targetBuffer.size());
  // TODO: we might want to find an appropriate name for this
  rearrangeData(sourceVector, data_dim, sourceBuffer);
  interpolator->interpolate(sourceData, targetData, data_dim);
  // TODO: same here
  rearrangeData(targetBuffer, data_dim, targetVector);
}

void read_displacement()
{
    // load the source mesh
    AMP::pout<<"Loading the source mesh"<<std::endl;
    boost::shared_ptr<AMP::InputDatabase> sourceMeshDatabase(new AMP::InputDatabase("SourceMesh"));
    sourceMeshDatabase->putString("MeshName", "mesh");
    sourceMeshDatabase->putString("MeshType", "libMesh");
    sourceMeshDatabase->putString("FileName", "/Users/qdi/Desktop/read_displacement/layer2_exodus.e");
    sourceMeshDatabase->putInteger("dim", 3);
    sourceMeshDatabase->putDouble("x_offset", 0.0);
    sourceMeshDatabase->putDouble("y_offset", 0.0);
    sourceMeshDatabase->putDouble("z_offset", 0.0);
    boost::shared_ptr<AMP::Mesh::MeshParameters> sourceMeshParams(new AMP::Mesh::MeshParameters(sourceMeshDatabase));
    sourceMeshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
    AMP::Mesh::Mesh::shared_ptr sourceMesh = AMP::Mesh::Mesh::buildMesh(sourceMeshParams);
    AMP::pout<<"***Loading the source mesh"<<std::endl;

    // build source vectors
    AMP::pout<<"Building the source vector"<<std::endl;
    bool const split = true;
    int const ghostWidth = 0;
    int const dofsPerNode = 3;
    AMP::Discretization::DOFManager::shared_ptr sourceDofManager = AMP::Discretization::simpleDOFManager::create(sourceMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
  AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("data"));
  AMP::LinearAlgebra::Vector::shared_ptr sourceVector = AMP::LinearAlgebra::createVector(sourceDofManager, variable, split);

std::cout<<"HERE"<<std::endl;
    std::string filename("/Users/qdi/Desktop/read_displacement/displacement");
    std::ifstream fin(filename.c_str());
    std::size_t n;
    fin>>n;
    ++n; // TODO:  I don't like that
    std::cout<<"Number of nodes ="<<n;
    std::vector<double> x(n), y(n), z(n);
    int node;
    for (std::size_t i = 0; i < n; ++i) {
        fin>>node>>x[i]>>y[i]>>z[i];
        std::cout<<i<<"  "<<node<<"  "<<std::setprecision(15)<<x[i]<<"  "<<y[i]<<"  "<<z[i]<<"\n";
    } // end for i
    if (!fin.eof()) std::cerr<<"WARNING:  didn't reach EOF"<<std::endl;
    fin.close();
}


void testDTK(AMP::UnitTest *ut, std::string exeName) {
  read_displacement();
abort();
  // build source and target meshes
  AMP::pout<<"Loading the source mesh"<<std::endl;
  boost::shared_ptr<AMP::InputDatabase> sourceMeshDatabase(new AMP::InputDatabase("SourceMesh"));
  sourceMeshDatabase->putString("MeshName", "mesh");
  sourceMeshDatabase->putString("MeshType", "libMesh");
  sourceMeshDatabase->putString("FileName", "/Users/qdi/Meshes/coupling.e");
  sourceMeshDatabase->putInteger("dim", 3);
  sourceMeshDatabase->putDouble("x_offset", 0.0);
  sourceMeshDatabase->putDouble("y_offset", 0.0);
  sourceMeshDatabase->putDouble("z_offset", 0.0);
  boost::shared_ptr<AMP::Mesh::MeshParameters> sourceMeshParams(new AMP::Mesh::MeshParameters(sourceMeshDatabase));
  sourceMeshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr sourceMesh = AMP::Mesh::Mesh::buildMesh(sourceMeshParams);
  
  AMP::pout<<"Building the target mesh"<<std::endl;
  boost::shared_ptr<AMP::InputDatabase> targetMeshDatabase(new AMP::InputDatabase("TargetMesh"));
  targetMeshDatabase->putString("MeshName", "mesh");
  targetMeshDatabase->putString("MeshType", "AMP");
  targetMeshDatabase->putInteger("dim", 3);
  targetMeshDatabase->putDouble("x_offset", 0.0);
  targetMeshDatabase->putDouble("y_offset", 0.0);
  targetMeshDatabase->putDouble("z_offset", 0.0);
  targetMeshDatabase->putString("Generator", "cube");
  std::vector<int> const size(3, 15);
  targetMeshDatabase->putIntegerArray("Size", size);
  std::vector<double> const range = sourceMesh->getBoundingBox();
  targetMeshDatabase->putDoubleArray("Range", range);
  boost::shared_ptr<AMP::Mesh::MeshParameters> targetMeshParams(new AMP::Mesh::MeshParameters(targetMeshDatabase));
  targetMeshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr targetMesh = AMP::Mesh::Mesh::buildMesh(targetMeshParams);

  // build dof managers
  bool const split = true;
  int const ghostWidth = 0;
  int const dofsPerNode = 3;
  AMP::Discretization::DOFManager::shared_ptr sourceDofManager = AMP::Discretization::simpleDOFManager::create(sourceMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);
  AMP::Discretization::DOFManager::shared_ptr targetDofManager = AMP::Discretization::simpleDOFManager::create(targetMesh, AMP::Mesh::Vertex, ghostWidth, dofsPerNode);

  // build source and target vectors
  AMP::LinearAlgebra::Variable::shared_ptr variable(new AMP::LinearAlgebra::Variable("data"));
  AMP::LinearAlgebra::Vector::shared_ptr sourceVector = AMP::LinearAlgebra::createVector(sourceDofManager, variable, split);
  AMP::LinearAlgebra::Vector::shared_ptr targetVector = AMP::LinearAlgebra::createVector(targetDofManager, variable, split);
  getVerticesCoordinates(sourceMesh->getIterator(AMP::Mesh::Vertex, ghostWidth), sourceDofManager, sourceVector);

  // build the interpolator and use it
  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm> const> commWrapper = Teuchos::opaqueWrapper<MPI_Comm>((sourceMesh->getComm()).getCommunicator());
  Teuchos::RCP<Teuchos::MpiComm<int> const> comm = Teuchos::createMpiComm<int>(commWrapper);
  // TODO: need to actually compute h
//double const  h = (range[1] - range[0]) / static_cast<double>(size[0]);
  double const h = 5.0;
  int const space_dim = 3;
  AMP::pout<<"source mesh size h = "<<h<<"\n";
  double const radius = 2.0*h;
//  std::string const interpolation_type = "Moving Least Square";
  std::string const interpolation_type = "Spline";
  std::string const basis_type = "Wendland";
  int const basis_order = 2;
  Teuchos::RCP<DataTransferKit::MeshFreeInterpolator> interpolator =
      DataTransferKit::MeshFreeInterpolatorFactory::create<int>(comm, interpolation_type, basis_type, basis_order, space_dim);

  AMP::pout<<"Setting up interpolator"<<std::endl;
  setupDTKInterpolator(sourceMesh, targetMesh, interpolator, radius/*should we pass it as an argument or infer it while traversing the source mesh?*/);
  
  AMP::pout<<"Performing the interpolation"<<std::endl;
  applyDTKInterpolator(interpolator, sourceVector, targetVector, 3);
  
  // postprocess
  targetVector->subtract(sourceVector, targetVector);
  AMP::pout<<"error l2 norm = "<<targetVector->L2Norm()<<"\n";
  
#ifdef USE_EXT_SILO
  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("silo");
  siloWriter->setDecomposition(1);
  siloWriter->registerVector(sourceVector, sourceMesh, AMP::Mesh::Vertex, "vector");
  siloWriter->writeFile("noname", 0);
  siloWriter.reset();
  siloWriter = AMP::Utilities::Writer::buildWriter("silo");
  siloWriter->setDecomposition(1);
  siloWriter->registerVector(targetVector, targetMesh, AMP::Mesh::Vertex, "vector");
  siloWriter->writeFile("noname", 1);
#endif

  ut->passes(exeName);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string exeName = "testDTK";

  try {
    testDTK(&ut, exeName);
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
