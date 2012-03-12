#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>
#include <limits>
#include <cmath>

#include "boost/shared_ptr.hpp"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"
#include "utils/ManufacturedSolution.h"

#include "ampmesh/SiloIO.h"

#include "ampmesh/Mesh.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/MultiVariable.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"

#include "libmesh.h"

#include "operators/diffusion/DiffusionNonlinearElement.h"
#include "operators/diffusion/DiffusionLinearElement.h"
#include "operators/diffusion/DiffusionLinearFEOperator.h"
#include "operators/diffusion/DiffusionNonlinearFEOperator.h"

#include "operators/boundary/NeumannVectorCorrection.h"

#include "../BVPOperatorParameters.h"
#include "../LinearBVPOperator.h"
#include "../NonlinearBVPOperator.h"
#include "../MassLinearFEOperator.h"
#include "../MassDensityModel.h"
#include "../OperatorBuilder.h"

#include "applyTests.h"


void bvpTest1(AMP::UnitTest *ut, const std::string exeName, const std::string meshName)
{
  // Tests diffusion Dirchlet BVP operator for temperature

  // Initialization
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logOnlyNodeZero(log_file);

  // Input database
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

//   Create the Mesh.
//--------------------------------------------------
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase(meshName.c_str());
  boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  boost::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh(mgrParams);
//--------------------------------------------------

  // Create nonlinear diffusion BVP operator and access volume nonlinear Diffusion operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> nonlinearPhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "FickNonlinearBVPOperator",
						   input_db,
						   nonlinearPhysicsModel);
  boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
          boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(nlinBVPOperator);
  boost::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nlinOp =
         boost::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(nlinBVPOp->getVolumeOperator());
  
  // use the linear BVP operator to create a linear diffusion operator with bc's
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> linearPhysicsModel;

  // Get source mass operator
  boost::shared_ptr<AMP::Operator::ElementPhysicsModel> sourcePhysicsModel;
  boost::shared_ptr<AMP::Operator::Operator> sourceOperator =
    AMP::Operator::OperatorBuilder::createOperator(meshAdapter,
						   "ManufacturedSourceOperator",
						   input_db,
						   sourcePhysicsModel);
  boost::shared_ptr<AMP::Operator::MassLinearFEOperator> sourceOp =
         boost::dynamic_pointer_cast<AMP::Operator::MassLinearFEOperator>(sourceOperator);


  boost::shared_ptr<AMP::Operator::MassDensityModel> densityModel = sourceOp->getDensityModel();
  boost::shared_ptr<AMP::ManufacturedSolution> mfgSolution = densityModel->getManufacturedSolution();

  // Set up input and output variables
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> tmp = 
      boost::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVariable>( nlinOp->getInputVariable() );
  boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> solVar(new AMP::LinearAlgebra::MultiVariable(tmp->getName()));
  for (size_t i=0; i<tmp->numVariables(); i++) {
      if ( tmp->getVariable(i).get() != NULL )
          solVar->add( tmp->getVariable(i) );
  }
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
  //boost::shared_ptr<AMP::Database> source_db = input_db->getDatabase("ManufacturedSourceOperator");
  //std::string sourceModelName = source_db->getString("LocalModel");
  //boost::shared_ptr<AMP::Database> sourceModel_db = input_db->getDatabase(sourceModelName);
  //boost::shared_ptr<AMP::Database> mfgSolution_db = sourceModel_db->getDatabase("ManufacturedSolution");
  /*bool isCylindrical = false;
  if (mfgSolution_db->keyExists("Geometry")) {
	  std::string geom = mfgSolution_db->getString("Geometry");
	  size_t pos = geom.find("Cylindrical");
	  size_t len = geom.size();
	  isCylindrical = (pos < len);
  }*/
  int zeroGhostWidth = 1;
  AMP::Mesh::MeshIterator  iterator = meshAdapter->getIterator(AMP::Mesh::Vertex, zeroGhostWidth);
  std::string mfgName = mfgSolution->get_name();
  if (mfgName.find("Cylindrical") < mfgName.size()) {
    for( ; iterator != iterator.end(); iterator++) {
		double x, y, z, r, th=0.;
    std::valarray<double> poly(10);
    x = ( iterator->coord() )[0];
    y = ( iterator->coord() )[1];
    z = ( iterator->coord() )[2];
		r = sqrt(x*x+y*y);
		double Pi=3.1415926535898;
		if (r>0) {th = acos(x/r); if (y<0.) th = 2*Pi-th;}
		mfgSolution->evaluate(poly,r,th,z);
    std::vector<size_t> gid;
    nodalDofMap->getDOFs ( iterator->globalID() , gid);
		solVec->setValueByGlobalID(gid[0], poly[0]);
	  }
  } else {
    for( ; iterator != iterator.end(); iterator++) {
		double x, y, z;
    std::valarray<double> poly(10);
    x = ( iterator->coord() )[0];
    y = ( iterator->coord() )[1];
    z = ( iterator->coord() )[2];
		mfgSolution->evaluate(poly,x,y,z);
    std::vector<size_t> gid;
    nodalDofMap->getDOFs ( iterator->globalID() , gid);
		solVec->setValueByGlobalID(gid[0], poly[0]);
	  }
  }

  solVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

  // Evaluate manufactured solution as an FE source
  sourceOp->apply(rhsVec, solVec, sourceVec, 1., 0.);

  // Evaluate action of diffusion operator
  nlinBVPOp->apply(sourceVec, solVec, resVec, 1., -1.);

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

      iterator = iterator.begin();
      size_t numNodes = 0;
      for(; iterator != iterator.end(); iterator++ ) numNodes++;

      iterator = iterator.begin();
      size_t iNode=0;
      double l2err = 0.;
      for(; iterator != iterator.end(); iterator++ ) {
        double x, y, z;
        x = ( iterator->coord() )[0];
        y = ( iterator->coord() )[1];
        z = ( iterator->coord() )[2];
        std::vector<size_t> gid;
        nodalDofMap->getDOFs ( iterator->globalID() , gid);
        double val, res, sol, src, err;
        res = resVec->getValueByGlobalID(gid[0]);
        sol = solVec->getValueByGlobalID(gid[0]);
        src = sourceVec->getValueByGlobalID(gid[0]);
        err = res/(src+.5*res + std::numeric_limits<double>::epsilon());
        std::valarray<double> poly(10);
        if (mfgName.find("Cylindrical") < mfgName.size()) {
    		double r = sqrt(x*x+y*y), th=0.;
    		double Pi=3.1415926535898;
    		if (r>0) {th = acos(x/r); if (y<0.) th = 2.*Pi-th;}
        	mfgSolution->evaluate(poly,r,th,z);
        } else {
        	mfgSolution->evaluate(poly,x,y,z);
        }
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
 #ifdef USE_SILO
     AMP::Mesh::SiloIO::shared_ptr  siloWriter( new AMP::Mesh::SiloIO);
     siloWriter->registerMesh( meshAdapter );

     siloWriter->registerVector(  workVec, meshAdapter, AMP::Mesh::Vertex, "RelativeError" );
     siloWriter->registerVector(   solVec, meshAdapter, AMP::Mesh::Vertex, "Solution" );
     siloWriter->registerVector(sourceVec, meshAdapter, AMP::Mesh::Vertex, "Source" );
     siloWriter->registerVector(   resVec, meshAdapter, AMP::Mesh::Vertex, "Residual" );
 
     siloWriter->writeFile( input_file , 0 );
 #endif

  ut->passes(exeName);
  std::cout.flush();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    std::vector<std::string> files, meshes;
    files.push_back("TensorDiffusion-Fick-MMS-1"); meshes.push_back("Mesh");
    files.push_back("TensorDiffusion-Fick-MMS-2"); meshes.push_back("Mesh");
    //files.push_back("Diffusion-Fick-OxMSRZC09-MMS-1");

    try {
        for (size_t i=0; i<files.size(); i++)
            bvpTest1(&ut, files[i], meshes[i]);
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


