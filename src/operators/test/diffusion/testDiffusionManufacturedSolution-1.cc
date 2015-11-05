#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <limits>
#include <cmath>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "utils/ManufacturedSolution.h"

#include "utils/Writer.h"

#include "ampmesh/Mesh.h"
#include "vectors/VectorBuilder.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "libmesh/libmesh.h"

#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/libmesh/NeumannVectorCorrection.h"

#include "operators/BVPOperatorParameters.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/libmesh/MassLinearFEOperator.h"
#include "operators/libmesh/MassDensityModel.h"
#include "operators/OperatorBuilder.h"

#include "../applyTests.h"


void bvpTest1(AMP::UnitTest *ut, std::string exeName)
{
  // Tests diffusion Dirchlet BVP operator for temperature

  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  // Input database
  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//--------------------------------------------------
//   Create the Mesh.
//--------------------------------------------------
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  AMP::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

  // Create nonlinear diffusion BVP operator and access volume nonlinear Diffusion operator
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> nonlinearPhysicsModel;
  AMP::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "ThermalNonlinearBVPOperator",
						   input_db,
						   nonlinearPhysicsModel);
  AMP::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          AMP::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nlinOp =
         AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  
  // use the linear BVP operator to create a linear diffusion operator with bc's
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> linearPhysicsModel;
  //AMP::shared_ptr<AMP::InputDatabase> bvp_db =
  //        AMP::dynamic_pointer_cast<AMP::InputDatabase>( input_db->getDatabase("ThermalLinearBVPOperator"));
  //AMP::shared_ptr<AMP::Operator::Operator> linBVPOperator =
  //        AMP::Operator::OperatorBuilder::createOperator(meshAdapter, bvp_db, linearPhysicsModel);
  //AMP::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
  //        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(linBVPOperator);
  //AMP::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> linOp =
  //        AMP::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(linBVPOp->getVolumeOperator());

  // Get source mass operator
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> sourcePhysicsModel;
  AMP::shared_ptr<AMP::Operator::Operator> sourceOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "ManufacturedSourceOperator",
						   input_db,
						   sourcePhysicsModel);
  AMP::shared_ptr<AMP::Operator::MassLinearFEOperator> sourceOp =
         AMP::dynamic_pointer_cast<AMP::Operator::MassLinearFEOperator>(sourceOperator);

  AMP::shared_ptr<AMP::Operator::MassDensityModel> densityModel = sourceOp->getDensityModel();
  AMP::shared_ptr<AMP::ManufacturedSolution> mfgSolution = densityModel->getManufacturedSolution();

  // Set up input and output vectors
  //AMP::LinearAlgebra::Variable::shared_ptr solVar = nlinOp->getInputVariable(nlinOp->getPrincipalVariableId());
  AMP::LinearAlgebra::Variable::shared_ptr solVar = nlinOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr rhsVar = nlinOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr resVar = nlinOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr sourceVar = sourceOp->getOutputVariable();
  AMP::LinearAlgebra::Variable::shared_ptr workVar = sourceOp->getOutputVariable();

  //----------------------------------------------------------------------------------------------------------------------------------------------//
  // Create a DOF manager for a nodal vector 
  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap = AMP::Discretization::simpleDOFManager::create(meshAdapter, AMP::Mesh::Vertex, nodalGhostWidth, DOFsPerNode, split);
  //----------------------------------------------------------------------------------------------------------------------------------------------//

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( nodalDofMap, solVar );
  AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( nodalDofMap, rhsVar );
  AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( nodalDofMap, resVar );
  AMP::LinearAlgebra::Vector::shared_ptr sourceVec = AMP::LinearAlgebra::createVector( nodalDofMap, sourceVar );
  AMP::LinearAlgebra::Vector::shared_ptr workVec = AMP::LinearAlgebra::createVector( nodalDofMap, workVar );

  rhsVec->setToScalar(0.0);

  // Fill in manufactured solution
  int zeroGhostWidth = 0;
  AMP::Mesh::MeshIterator  iterator = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
  for( ; iterator != iterator.end(); ++iterator)
  {
    double x, y, z;
    std::valarray<double> poly(10);
    x = ( iterator->coord() )[0];
    y = ( iterator->coord() )[1];
    z = ( iterator->coord() )[2];
    mfgSolution->evaluate(poly,x,y,z);
    std::vector<size_t> i;
    nodalDofMap->getDOFs ( iterator->globalID() , i);
    solVec->setValueByGlobalID(i[0], poly[0]);
  }
  solVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  // Evaluate manufactured solution as an FE source
  sourceOp->apply(solVec, sourceVec);

  // Evaluate action of diffusion operator
  nlinBVPOp->residual(sourceVec, solVec, resVec);
  resVec->scale(-1.0);

  // Output Mathematica form (requires serial execution)
  for (int i=0; i<globalComm.getSize(); i++) {
    if ( globalComm.getRank()==i ) {
      std::string filename="data_"+exeName;
      int rank = globalComm.getRank();
      int nranks = globalComm.getSize();
      std::ios_base::openmode omode=std::ios_base::out;
      if (rank>0) omode |= std::ios_base::app;
      std::ofstream file(filename.c_str(),omode);
      if (rank == 0) {
          file << "(* x y z solution solution fe-source fe-operator error *)" << std::endl;
          file << "results={" << std::endl;
      }

      AMP::Mesh::MeshIterator  urIterator = iterator.begin();
      size_t numNodes = 0;
      for(; urIterator != iterator.end(); urIterator++ ) numNodes++;

      size_t iNode=0;
      double l2err = 0.;
      AMP::Mesh::MeshIterator  myIterator = iterator.begin();
      for(; myIterator != iterator.end(); myIterator++ ) {
        double x, y, z;
        x = ( myIterator->coord() )[0];
        y = ( myIterator->coord() )[1];
        z = ( myIterator->coord() )[2];
        std::vector<size_t> gid;
        nodalDofMap->getDOFs ( myIterator->globalID() , gid);
        double val, res, sol, src, err;
        res = resVec->getValueByGlobalID(gid[0]);
        sol = solVec->getValueByGlobalID(gid[0]);
        src = sourceVec->getValueByGlobalID(gid[0]);
        err = res/(src+.5*res + std::numeric_limits<double>::epsilon());
        std::valarray<double> poly(10);
        mfgSolution->evaluate(poly,x,y,z);
        val = poly[0];
        workVec->setValueByGlobalID(gid[0], err);

        file << "{" << x << "," << y << "," << z <<"," << val <<  ","
                << sol << "," << src << "," << res+src << "," << err << "}";
        if (iNode<numNodes-1) file << "," << std::endl;

        l2err += (res*res);
        iNode++;
      }

      if (rank == nranks-1) {
          file << "};" << std::endl;
          file << "nodes = " << numNodes <<"; l2err = " << l2err << ";" << std::endl;
      }

      file.close();
    }
    globalComm.barrier();
  }

  // Plot the results
  if( globalComm.getSize() == 1 ) {
 #ifdef USE_EXT_SILO
     AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector(  workVec, meshAdapter, AMP::Mesh::Vertex, "RelativeError" );
     siloWriter->registerVector(   solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
     siloWriter->registerVector(sourceVec, meshAdapter, AMP::Mesh::Vertex, "Source" );
     siloWriter->registerVector(   resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
 
     siloWriter->writeFile( input_file , 0 );
 #endif
   }

  ut->passes(exeName);
  std::cout.flush();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    std::vector<std::string> files;
    files.push_back("Diffusion-Fick-TUI-MMS-1");
    //files.push_back("Diffusion-Fick-OxMSRZC09-MMS-1");

    try {
        for (size_t i=0; i<files.size(); i++)
            bvpTest1(&ut, files[i]);
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


