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

const size_t dofs_per_var = 10; // dofs per variable; number of axial faces
const size_t num_dofs = 2*dofs_per_var; // total number of dofs

// function to check that Jacobian matches known values
bool JacobianIsCorrect(boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian, double knownJacobian[num_dofs][num_dofs])
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

  // set dummy values for reset in operator constructor; otherwise zero-values give error in thermodynamic property evaluations
  for (; face != face.end(); face++){
     faceDOFManager->getDOFs( face->globalID(), dofs );
     FrozenVec->setValueByGlobalID(dofs[0], 900.0e3);
     FrozenVec->setValueByGlobalID(dofs[1], 15.0e6);;
     SolVec->setValueByGlobalID(dofs[0],1.0);
     SolVec->setValueByGlobalID(dofs[1],1.0);
  }

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

  // reset the linear operator
  subchannelOperator->reset(subchannelOpParams);

   {//test block #1
      std::cout<<std::endl<<"Test #1 Computed Jacobian:"<<std::endl;
      std::vector<size_t> dofs;
      AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],1000.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 16.4e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],900.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 16.3e6);;
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
      FrozenVec->setValueByGlobalID(dofs[0],700.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 16.1e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],300.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 13.5e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],450.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 9.0e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],570.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 12.0e5);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],230.0e2);
      FrozenVec->setValueByGlobalID(dofs[1], 4.0e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],999.9e3);
      FrozenVec->setValueByGlobalID(dofs[1], 14.0e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],235.6e3);
      FrozenVec->setValueByGlobalID(dofs[1], 12.5e6);
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
    
      subchannelOperator->setFrozenVector(FrozenVec);
      subchannelOperator->reset(subchannelOpParams);
      subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    
      // get the matrix
      boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();
    
      double knownJacobian[num_dofs][num_dofs] = {
         {1,0.000939452905154014,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0.00124759679718839,-1.00000429530722,0.011268692736763,0.999960685635811,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0.00105568179062782,-1.00000368307659,0.00999584631218802,0.999964874076558,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0.000913511255920129,-1.00000321012603,0.00899327482841767,0.999968798175212,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0.000776481818973986,-1.00000269397412,0.00462823638551236,0.999981493118125,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0.000388303980414046,-1.00000155270719,0.00669650537073264,0.999977181456216,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0.000585817765508914,-1.00000199619168,0.00851625738946539,0.999969145416505,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0.000721218570384414,-1.00000261299038,0.00125280086469404,0.999985692192925,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000113887289793515,-1.0000013006675,0.0130211335246558,0.999954800323239,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00118338425166691,-1.00000410782867,0.00379642355524263,0.999982176930523},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}
      };
      bool passedJacobianTest = JacobianIsCorrect(testJacobian,knownJacobian);
      if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #1");
      else ut->failure(exeName+": apply: known Jacobian value test #1");
   }//end of test block #1

   {//test block #2
      std::cout<<std::endl<<"Test #2 Computed Jacobian:"<<std::endl;
      std::vector<size_t> dofs;
      AMP::Mesh::MeshIterator face = xyFaceMesh->getIterator(AMP::Mesh::Face, 0);
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],950.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 15.0e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],850.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 15.1e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],700.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 15.25e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],500.0e3);
      FrozenVec->setValueByGlobalID(dofs[1], 15.26e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],324.6e3);
      FrozenVec->setValueByGlobalID(dofs[1], 11.0e5);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],457.7e3);
      FrozenVec->setValueByGlobalID(dofs[1], 12.5e5);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],134.6e2);
      FrozenVec->setValueByGlobalID(dofs[1], 34.5e5);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],457.6e3);
      FrozenVec->setValueByGlobalID(dofs[1], 12.0e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],325.7e3);
      FrozenVec->setValueByGlobalID(dofs[1], 11.5e6);;
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
      ++face;
      faceDOFManager->getDOFs( face->globalID(), dofs );
      FrozenVec->setValueByGlobalID(dofs[0],898.6e3);
      FrozenVec->setValueByGlobalID(dofs[1], 15.7e6);
      SolVec->setValueByGlobalID(dofs[0],1.0);
      SolVec->setValueByGlobalID(dofs[1],1.0);
    
      subchannelOperator->setFrozenVector(FrozenVec);
      subchannelOperator->reset(subchannelOpParams);
      subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    
      // get the matrix
      boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();
    
      double knownJacobian[num_dofs][num_dofs] = {
         {1,0.00094733031302362,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0.00115214013912312,-1.00000401532523,0.0106360680973219,0.999962385419614,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0.000978708909950683,-1.00000346121561,0.00903331568143602,0.999968413640661,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0.000795304332477446,-1.00000278090231,0.00708702135053376,0.999976234623364,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0.00060005764277089,-1.00000201221291,0.00520094270260578,0.999980754106999,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0.000442342995882102,-1.00000163687363,0.00716016681591477,0.999975697667819,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0.00059706187256616,-1.00000202648853,0.00118533412239022,0.999986316566661,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,9.84847789253265e-05,-1.00000113690299,0.00666789816034634,0.99997734081775,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000563297800467833,-1.00000191422652,0.00499264551266783,0.999981035371265,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000450602583900622,-1.00000171161976,0.0111938850587298,0.9999608026051},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}
      };
      bool passedJacobianTest = JacobianIsCorrect(testJacobian,knownJacobian);
      if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #2");
      else ut->failure(exeName+": apply: known Jacobian value test #2");
   }//end of test block #2

   {//test block #3
      std::cout<<std::endl<<"Test #3 Computed Jacobian:"<<std::endl;
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
      subchannelOperator->reset(subchannelOpParams);
      subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);
    
      // get the matrix
      boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();
    
      double knownJacobian[num_dofs][num_dofs] = {
         {1,0.000962166060627425,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0.000853056733520424,-1.00000305938232,0.0114067382625288,0.999959276096689,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0.00107453717325675,-1.00000383627177,0.00999779153622338,0.99996486724094,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0.00091397270775785,-1.0000032117476,0.00919870562562558,0.999964753297376,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0.000802428956144871,-1.00000307466897,0.00580374927940146,0.999979600978116,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0.00050585549305329,-1.00000177798123,0.00871725609681695,0.999969151111928,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0.000799562175628259,-1.00000282951468,0.00997834397602954,0.999963938940236,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0.000863464422657599,-1.00000312050198,0.000943493036102408,0.999983499948451,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,8.03714696059628e-05,-1.00000140555716,0.00879232188379425,0.999969023095884,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.000767452302123088,-1.00000270387011,0.00567135482253266,0.999979886341977},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,1,0},
         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}
      };
      bool passedJacobianTest = JacobianIsCorrect(testJacobian,knownJacobian);
      if (passedJacobianTest) ut->passes(exeName+": apply: known Jacobian value test #3");
      else ut->failure(exeName+": apply: known Jacobian value test #3");
   }//end of test block #3
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


