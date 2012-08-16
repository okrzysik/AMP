#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

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
#include "operators/OperatorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"

const size_t num_cells_per_var = 9;
const size_t num_faces_per_var = num_cells_per_var+1;
const size_t num_faces = 2*num_faces_per_var;

//
// function to check that Jacobian matches known values
bool JacobianIsCorrect(boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian, double knownJacobian[num_faces][num_faces])
{
   bool passed = true; // boolean for all values being equal to known values
   //double max_abs = 0.0; // maximum absolute error; L-infinity norm of error
   //double max_rel = 0.0; // maximum relative error
   
/*
   int DofsPerFace =  2;
   AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( xyFaceMesh, AMP::Mesh::Face, 1, DofsPerFace, true);
   AMP::Mesh::MeshIterator iface = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
   AMP::Mesh::MeshIterator jface = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
   std::vector<size_t> idofs;
   std::vector<size_t> jdofs;
   for (size_t i = 0; i <= 9; ++iface, ++i){
      faceDOFManager->getDOFs( iface->globalID(), idofs );
      for (size_t j = 0; j <= 9; ++jface, ++j){
         faceDOFManager->getDOFs( jface->globalID(), jdofs );
         // get Jacobian values
         double J00 = testJacobian->getValueByGlobalID(idofs[0],jdofs[0]);
         double J01 = testJacobian->getValueByGlobalID(idofs[0],jdofs[1]);
         double J10 = testJacobian->getValueByGlobalID(idofs[1],jdofs[0]);
         double J11 = testJacobian->getValueByGlobalID(idofs[1],jdofs[1]);
         // check that Jacobian value is approximately equal to known value
         double known00 = knownJacobian[2*i  ][2*j  ];
         double known01 = knownJacobian[2*i  ][2*j+1];
         double known10 = knownJacobian[2*i+1][2*j  ];
         double known11 = knownJacobian[2*i+1][2*j+1];
         if (!AMP::Utilities::approx_equal(J00,known00,0.01)) passed = false;
         if (!AMP::Utilities::approx_equal(J01,known01,0.01)) passed = false;
         if (!AMP::Utilities::approx_equal(J10,known10,0.01)) passed = false;
         if (!AMP::Utilities::approx_equal(J11,known11,0.01)) passed = false;
         // calculate absolute error
         double abs00 = std::abs(J00 - known00);
         double abs01 = std::abs(J01 - known01);
         double abs10 = std::abs(J10 - known10);
         double abs11 = std::abs(J11 - known11);
         max_abs = std::max(max_abs,abs00);
         max_abs = std::max(max_abs,abs01);
         max_abs = std::max(max_abs,abs10);
         max_abs = std::max(max_abs,abs11);
         // calculate relative error
         double rel00;
         double rel01;
         double rel10;
         double rel11;
         if (std::abs(known00) < 1.0e-15) rel00 = 0.0;
         else rel00 = abs00 / known00;
         if (std::abs(known01) < 1.0e-15) rel01 = 0.0;
         else rel01 = abs01 / known01;
         if (std::abs(known10) < 1.0e-15) rel10 = 0.0;
         else rel10 = abs10 / known10;
         if (std::abs(known11) < 1.0e-15) rel11 = 0.0;
         else rel11 = abs11 / known11;
         max_rel = std::max(max_rel,rel00);
         max_rel = std::max(max_rel,rel01);
         max_rel = std::max(max_rel,rel10);
         max_rel = std::max(max_rel,rel11);
      }
   }
*/
      for(size_t i = 0; i < num_faces; i++) {
          std::vector<unsigned int> matCols;
          std::vector<double> matVals;
          testJacobian->getRowByGlobalID(i, matCols, matVals);
          std::cout<<"{";
          size_t col=0;
          // for nonzero entries of row i
          for(size_t j = 0; j < matCols.size(); j++) {
            // zeros before first nonzero
            while (col<matCols[j]) {
                std::cout << "0";
                std::cout<<",";
                if (!AMP::Utilities::approx_equal(0.0,knownJacobian[i][col],0.01)) passed = false;
                col++;
            }
            std::cout<<matVals[j];
            if (!AMP::Utilities::approx_equal(matVals[j],knownJacobian[i][col],0.01)) passed = false;
            if (matCols[j]<num_faces-1) std::cout<<",";
            col++;
          }//end for j
          // zeros after last nonzero
          while (col<num_faces) {
            std::cout << "0";
            if (col<num_faces-1)std::cout<<",";
            if (!AMP::Utilities::approx_equal(0.0,knownJacobian[i][col],0.01)) passed = false;
            col++;
          }
          std::cout<<"}";
          if (i<num_faces-1) std::cout<<","<<std::endl;
      }//end for i

/*
   if (!passed){
      // report errors
      AMP::pout << "Jacobian did not match known Jacobian.\n";
      AMP::pout << "Jacobian comparison: Maximum absolute error: " << max_abs << ", Maximum relative error: " << max_rel << "\n";

      AMP::pout << "Test Jacobian:\n";
      for (size_t i = 0; i <= 9; i++){
         AMP::pout << "\n";
         for (size_t j = 0; j <= 9; j++){
            AMP::pout << "   " << testJacobian[i][j];
         }
      }
      AMP::pout << "\n\nKnown Jacobian:\n";
      for (size_t i = 0; i <= 9; i++){
         AMP::pout << "\n";
         for (size_t j = 0; j <= 9; j++){
            AMP::pout << "   " << knownJacobian[i][j];
         }
      }
   }
*/
   return passed;
}

void Test(AMP::UnitTest *ut, const std::string exeName)
{
  // create input and output file names
  std::string input_file = "input_"  + exeName;
  std::string log_file   = "output_" + exeName;
  AMP::PIO::logOnlyNodeZero(log_file);

  // get input database from input file
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // create mesh
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> meshParams(new AMP::Mesh::MeshParameters(mesh_db));
  meshParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  boost::shared_ptr<AMP::Mesh::Mesh> subchannelMesh = AMP::Mesh::Mesh::buildMesh(meshParams);
  AMP::Mesh::Mesh::shared_ptr xyFaceMesh;
  xyFaceMesh = subchannelMesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator( subchannelMesh , 0 ) );

  // get dof manager
  int DofsPerFace =  2;
  AMP::Discretization::DOFManager::shared_ptr faceDOFManager = AMP::Discretization::simpleDOFManager::create( subchannelMesh, 
		  AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1), AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0), DofsPerFace);

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
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],700.0e3);
  FrozenVec->setValueByGlobalID(dofs[1], 12.4e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],900.0e3);
  FrozenVec->setValueByGlobalID(dofs[1], 12.3e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],800.0e3);
  FrozenVec->setValueByGlobalID(dofs[1], 16.2e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],650.0e3);
  FrozenVec->setValueByGlobalID(dofs[1], 14.1e5);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],367.4e3);
  FrozenVec->setValueByGlobalID(dofs[1], 31.5e5);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],657.2e3);
  FrozenVec->setValueByGlobalID(dofs[1], 12.5e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],788.5e3);
  FrozenVec->setValueByGlobalID(dofs[1], 12.7e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],235.7e2);
  FrozenVec->setValueByGlobalID(dofs[1], 17.8e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],673.1e3);
  FrozenVec->setValueByGlobalID(dofs[1], 13.6e6);;
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);
  ++face;
  faceDOFManager->getDOFs( face->globalID(), dofs );
  FrozenVec->setValueByGlobalID(dofs[0],385.2e3);
  FrozenVec->setValueByGlobalID(dofs[1], 16.3e6);
  SolVec->setValueByGlobalID(dofs[0],1.0);
  SolVec->setValueByGlobalID(dofs[1],1.0);

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // create linear operator
  // get linear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelTwoEqLinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = subchannelMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
  subchannelOpParams->d_frozenSolution = FrozenVec ;
  subchannelOpParams->d_dofMap = faceDOFManager ;

  boost::shared_ptr<AMP::Operator::SubchannelTwoEqLinearOperator> subchannelOperator (new AMP::Operator::SubchannelTwoEqLinearOperator(subchannelOpParams));


  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);

  /*
     {
  // Test apply with known residual evaluation
  FrozenVec->setValueByLocalID(0,1000.0);
  FrozenVec->setValueByLocalID(1,15.3e6);
  FrozenVec->setValueByLocalID(2,15.2e6);
  FrozenVec->setValueByLocalID(3,15.1e6);
  FrozenVec->setValueByLocalID(4,15.0e6);
  SolVec->setValueByLocalID(0,1.0);
  SolVec->setValueByLocalID(1,1.0);
  SolVec->setValueByLocalID(2,1.0);
  SolVec->setValueByLocalID(3,1.0);
  SolVec->setValueByLocalID(4,1.0);
  subchannelOperator->setFrozenVector(FrozenVec);
  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

  double knownJacobian[20][20] = {
  {1.00000286783759,0.000945636121207302,0,0,0},
  {0.0084145730631915,-1.00006033273598,0.999924832234253,0,0},
  {0.0164760871167386,0,-1.00006231804521,0.999922391414922,0},
  {0.0240080126558191,0,0,-1.00006376793541,0.999920554193017},
  {0,0,0,0,0}
  };
  bool passedJacobianTest ;
  //     bool passedJacobianTest = JacobianIsCorrect(subchannelOperator->d_Jacobian,knownJacobian);
  if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #1");
  else ut->failure(exeName+": apply: known Jacobian value test #1");
  }

  {
  // Test apply with known residual evaluation
  FrozenVec->setValueByLocalID(0,500.0e3);
  FrozenVec->setValueByLocalID(1,17.0e6);
  FrozenVec->setValueByLocalID(2,17.1e6);
  FrozenVec->setValueByLocalID(3,17.2e6);
  FrozenVec->setValueByLocalID(4,17.3e6);
  SolVec->setValueByLocalID(0,1.0);
  SolVec->setValueByLocalID(1,1.0);
  SolVec->setValueByLocalID(2,1.0);
  SolVec->setValueByLocalID(3,1.0);
  SolVec->setValueByLocalID(4,1.0);
  subchannelOperator->setFrozenVector(FrozenVec);
  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

  double knownJacobian[5][5] = {
  {1.00000000036425,0.00093610492591568,0,0,0},
  {0.0616135012140977,-1.00009118544544,0.999886751731324,0,0},
  {0.068075038946247,0,-1.00010127215554,0.99987438016631,0},
  {0.073216467438169,0,0,-1.00011121418204,0.999862069911393},
  {0,0,0,0,0}
  };
  bool passedJacobianTest ; 
  //     bool passedJacobianTest = JacobianIsCorrect(subchannelOperator->d_Jacobian,knownJacobian);
  if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #2");
  else ut->failure(exeName+": apply: known Jacobian value test #2");
  }
  */

  {
	  std::vector<size_t> dofs;
	  AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],700.0e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 12.4e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],900.0e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 12.3e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],800.0e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 16.2e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],650.0e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 14.1e5);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],367.4e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 31.5e5);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],657.2e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 12.5e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],788.5e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 12.7e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],235.7e2);
	  FrozenVec->setValueByGlobalID(dofs[1], 17.8e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],673.1e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 13.6e6);;
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);
	  ++face;
	  faceDOFManager->getDOFs( face->globalID(), dofs );
	  FrozenVec->setValueByGlobalID(dofs[0],385.2e3);
	  FrozenVec->setValueByGlobalID(dofs[1], 16.3e6);
	  SolVec->setValueByGlobalID(dofs[0],1.0);
	  SolVec->setValueByGlobalID(dofs[1],1.0);

	  subchannelOperator->setFrozenVector(FrozenVec);
	  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

	  // get the matrix
	  boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();

/*
	  double knownJacobian[num_faces][num_faces] = {
		  {0.99999999749821,0.000962166060627425,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0.00138794996064453,-1.00000518247891,0.0107107344615788,0.999961470497116,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {-0.99999999749821,0,0.999999999996379,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0.000657640227042369,-1.00000227648592,0.0103580346127423,0.999963346640145,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,-0.999999999996379,0,1.00000000120429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0.00066102699712792,-1.00000336910008,0.00936730839615044,0.999964702837787,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,-1.00000000120429,0,0.999999998049371,0,0,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,5.20602361744438e-05,-0.999998655019239,0.00687742037746793,0.999976031436162,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,-0.999999998049371,0,1.00000000005676,0,0,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0.00144432982241215,-1.000004748495,0.00807042369432943,0.999971857074367,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,-1.00000000005676,0,1.00000000033991,0,0,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0.00111498405002515,-1.0000042409559,0.00962716112679294,0.999965292605162,0,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,-1.00000000033991,0,1.00000000215654,0,0,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,-0.00131362150614084,-0.999994820859473,0.00365092014412878,0.999981263374941,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,-1.00000000215654,0,0.999999928144391,0,0,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00238396838980072,-1.00000288918105,0.00683624602531728,0.99997541183633,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.999999928144391,0,1.00000000271785,0,0,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.25685271322929e-06,-0.999999346190506,0.00661590812865051,0.999977304559597},
		  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.00000000271785,0,1.00000000033377,0},
		  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.999999998997306}
	  };
*/
	  double knownJacobian[num_faces][num_faces] = {
      {0.99999999749821,0.000962166060627425,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0.000853057022888615,-1.00000306169866,0.0114067381600622,0.999959277908827,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {-0.99999999749821,0,0.999999999996379,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0.00107453511147884,-1.00000383762854,0.00999779404243742,0.999964866457383,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,-0.999999999996379,0,1.00000000120429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0.000914010820042994,-1.00000320911932,0.00919875828425407,0.99996478638668,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,-1.00000000120429,0,0.999999998049371,0,0,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0.000802425398302349,-1.00000307788877,0.00580375217670267,0.999979600616436,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,-0.999999998049371,0,1.00000000005676,0,0,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0.000505812031269878,-1.00000179404584,0.00871725140241169,0.999969147582058,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,-1.00000000005676,0,1.00000000033991,0,0,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0.00079956227573068,-1.00000282760225,0.0099783445489746,0.999963938105396,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,-1.00000000033991,0,1.00000000215654,0,0,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0.000863471765111286,-1.00000311815791,0.000944591248788966,0.999983501982161,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,-1.00000000215654,0,0.999999928144391,0,0,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,7.99653967228754e-05,-1.00000140670354,0.00879232195117482,0.999969021401661,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.999999928144391,0,1.00000000271785,0,0,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000767470804310478,-1.00000270273871,0.00567136820097811,0.999979885587762},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.00000000271785,0,1.00000000033377,0},
      {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.999999998997306}
          };
	  bool passedJacobianTest = JacobianIsCorrect(testJacobian,knownJacobian);
	  if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #3");
	  else ut->failure(exeName+": apply: known Jacobian value test #3");
  }

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


