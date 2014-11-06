#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "utils/shared_ptr.h"

#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/PIO.h"

#include "vectors/VectorBuilder.h"

#include "operators/subchannel/SubchannelPhysicsModel.h"
#include "operators/subchannel/SubchannelOperatorParameters.h"
#include "operators/subchannel/SubchannelTwoEqLinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/structuredFaceDOFManager.h"

const size_t dofs_per_var = 10; // dofs per variable; number of axial faces
const size_t num_dofs = 2*dofs_per_var; // total number of dofs

// function to check that Jacobian matches known values
bool JacobianIsCorrect(AMP::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian, double knownJacobian[num_dofs][num_dofs])
{
   bool passed = true; // boolean for all values being equal to known values
   std::stringstream mismatch; // string containing error messages for mismatched Jacobian entries
      
   // loop over rows of Jacobian
   for(size_t i = 0; i < num_dofs; i++) {
       std::vector<unsigned int> matCols; // indices of nonzero entries in row i
       std::vector<double> matVals;       // values of nonzero entries in row i
       testJacobian->getRowByGlobalID(i, matCols, matVals); // get nonzero entries of row i of Jacobian
       std::cout<<"{";
       size_t col=0;
       // loop over nonzero entries of row i
       for(size_t j = 0; j < matCols.size(); j++) {
         // zeros before first nonzero entry
         while (col<matCols[j]) {
             std::cout << "0";
             std::cout<<",";
             if (!AMP::Utilities::approx_equal(0.0,knownJacobian[i][col],0.01)){
                passed = false;
                mismatch<<"Entry does not match. i = "<<i<<", j = "<<col<<", Computed = 0.0, Known = "<<knownJacobian[i][col]<<std::endl;
             }
             col++;
         }
         std::cout<<matVals[j];
         if (!AMP::Utilities::approx_equal(matVals[j],knownJacobian[i][col],0.01)){
            passed = false;
            mismatch<<"Entry does not match. i = "<<i<<", j = "<<col<<", Computed = "<<matVals[j]<<", Known = "<<knownJacobian[i][col]<<std::endl;
         }
         if (matCols[j]<num_dofs-1) std::cout<<",";
         col++;
       }//end for j
       // zeros after last nonzero entry
       while (col<num_dofs) {
         std::cout << "0";
         if (col<num_dofs-1)std::cout<<",";
         if (!AMP::Utilities::approx_equal(0.0,knownJacobian[i][col],0.01)){
            passed = false;
            mismatch<<"Entry does not match. i = "<<i<<", j = "<<col<<", Computed = 0.0, Known = "<<knownJacobian[i][col]<<std::endl;
         }
         col++;
       }
       std::cout<<"}";
       if (i<num_dofs-1) std::cout<<","<<std::endl;
   }//end for i
   std::cout<<std::endl;
   std::cout<<mismatch.str();

   return passed;
}

void Test(AMP::UnitTest *ut, std::string exeName)
{
  // create input and output file names
  std::string input_file = "input_"  + exeName;
  std::string log_file   = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);

  // get input database from input file
  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // create mesh
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  AMP::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
  xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

  // get dof manager
  int DOFsPerFace[3]={0,0,2};
  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = 
    AMP::Discretization::structuredFaceDOFManager::create( subchannelMesh, DOFsPerFace, 1 );

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  (new AMP::LinearAlgebra::Variable("flow"));
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable  (new AMP::LinearAlgebra::Variable("flow"));

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr FrozenVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr SolVec = AMP::LinearAlgebra::createVector( faceDOFManager , inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr ResVec = AMP::LinearAlgebra::createVector( faceDOFManager , outputVariable, true );

  // set frozen vector before construction of linear operator to prevent reset from being applied
  // with a zero frozen vector
  std::vector<size_t> dofs;
  AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);

  const double h_scale = AMP::Operator::Subchannel::scaleEnthalpy;
  const double P_scale = AMP::Operator::Subchannel::scalePressure;

  // set dummy values for reset in operator constructor; otherwise zero-values give error in thermodynamic property evaluations
  for (; face != face.end(); face++){
     faceDOFManager->getDOFs( face->globalID(), dofs );
     FrozenVec->setValueByGlobalID(dofs[0], h_scale*900.0e3);
     FrozenVec->setValueByGlobalID(dofs[1], P_scale*15.0e6);;
     SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
     SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
  }

  // create subchannel physics model
  AMP::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  AMP::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( 
    new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  AMP::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel(
    new AMP::Operator::SubchannelPhysicsModel(params));

  // create linear operator
  // get linear operator database
  AMP::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelTwoEqLinearOperator");
  // set operator parameters
  AMP::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(
    new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = subchannelMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
  subchannelOpParams->d_frozenSolution = FrozenVec ;
  subchannelOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
  subchannelOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
  subchannelOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");
  AMP::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelOperator(
    new AMP::Operator::SubchannelTwoEqLinearOperator(subchannelOpParams));

  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // reset the linear operator
  subchannelOperator->reset(subchannelOpParams);

   {//test block
      std::cout<<std::endl<<"Test Computed Jacobian:"<<std::endl;
      std::vector<size_t> dofs;
      AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*700.0e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 12.4e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*900.0e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 12.3e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*800.0e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 16.2e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*650.0e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 14.1e5);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*367.4e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 31.5e5);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*657.2e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 12.5e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*788.5e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 12.7e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*235.7e2);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 17.8e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*673.1e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 13.6e6);;
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],h_scale*385.2e3);
      FrozenVec->setValueByGlobalID(dofs[1],P_scale* 16.3e6);
      SolVec->setValueByGlobalID(dofs[0],h_scale*1.0);
      SolVec->setValueByGlobalID(dofs[1],P_scale*1.0);
    
      subchannelOperator->setFrozenVector(FrozenVec);
      subchannelOpParams->d_initialize = true;
      subchannelOperator->reset(subchannelOpParams);
      subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    
      // get the matrix
      AMP::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();
    
      double knownJacobian[num_dofs][num_dofs] = {
         {0.99999999749821,0.000778407237773775,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {-3.35541872265853e-07,-8.78766222791986e-05,3.2990887838029e-07,8.78766284349374e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {-0.767213114757281,0,0.767213115083247,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,-3.51799478071037e-07,-8.78765679184489e-05,3.45890251361492e-07,8.78765921713374e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,-0.767213115083247,0,0.767213116091066,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,-3.12341330405279e-07,-8.78767302051646e-05,3.15095630760153e-07,8.78765973055314e-05,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,-0.767213115515851,0,0.767213114998425,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,-2.96651446897112e-07,-8.78767082179109e-05,1.91141714652964e-07,8.78771187267925e-05,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,-0.767213114998425,0,0.767213113877773,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,-1.91535569085153e-07,-8.7877166661357e-05,2.91551133926647e-07,8.78767591853259e-05,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,-0.767213116382785,0,0.767213115308738,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,-2.75318671573501e-07,-8.78768536463993e-05,3.41860441331691e-07,8.787656753596e-05,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,-0.767213116008937,0,0.767213114510477,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,-3.25732837135205e-07,-8.78766789006946e-05,1.85314096485698e-08,8.78772557339642e-05,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,-0.767213114510477,0,0.767213005345743,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.64429155390083e-08,-8.7877285625751e-05,2.86792534131204e-07,8.78767567742423e-05,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.767213083436951,0,0.767213116597182,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3.35566038867217e-07,-8.78766570500576e-05,1.54853676092429e-07,8.78772425610155e-05},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.767213116597182,0,0.767213114375646,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.999999998997306}
      };
      bool passedJacobianTest = JacobianIsCorrect(testJacobian,knownJacobian);
      if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test");
      else ut->failure(exeName+": apply: known Jacobian value test");
   }//end of test block
}

int main(int argc, char *argv[])
{
	AMP::AMPManagerProperties startup_properties;
	startup_properties.use_MPI_Abort = false;
	AMP::AMPManager::startup(argc,argv,startup_properties);

	AMP::UnitTest ut;

	const int NUMFILES=1;
	std::string files[NUMFILES] = {
		"testSubchannelTwoEqLinearOperator"
	};

	for (int i=0; i<NUMFILES; i++) {
		try {
			Test(&ut, files[i]);

		} catch (std::exception &err) {
			std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
			ut.failure("ERROR: While testing: "+files[i]);
		} catch( ... ) {
			std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
			ut.failure("ERROR: While testing: "+files[i]);
		}
	}

	ut.report();

	int num_failed = ut.NumFailGlobal();
	AMP::AMPManager::shutdown();
	return num_failed;
}   


