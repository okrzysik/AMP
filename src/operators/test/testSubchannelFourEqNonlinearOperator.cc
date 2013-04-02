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
#include "operators/subchannel/SubchannelFourEqNonlinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"

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
  AMP::Discretization::DOFManager::shared_ptr subchannelDOFManager;
  if ( subchannelMesh.get() != NULL ) {
    AMP::Mesh::MeshIterator axialFaces0 = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,0);
    AMP::Mesh::MeshIterator axialFaces1 = AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(subchannelMesh,1);
    AMP::Mesh::MeshIterator gapFaces0 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, 
      AMP::Mesh::StructuredMeshHelper::getXZFaceIterator(subchannelMesh,0), 
      AMP::Mesh::StructuredMeshHelper::getYZFaceIterator(subchannelMesh,0) );
    AMP::Mesh::MeshIterator gapFaces1 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, 
      AMP::Mesh::StructuredMeshHelper::getXZFaceIterator(subchannelMesh,1), 
      AMP::Mesh::StructuredMeshHelper::getYZFaceIterator(subchannelMesh,1) );
    std::vector<AMP::Discretization::DOFManager::shared_ptr> subchannelChildrenDOFManagers(2);
    subchannelChildrenDOFManagers[0] = AMP::Discretization::simpleDOFManager::create( subchannelMesh, axialFaces1, axialFaces0, 3 );
    subchannelChildrenDOFManagers[1] = AMP::Discretization::simpleDOFManager::create( subchannelMesh, gapFaces1, gapFaces0, 1 );
    subchannelDOFManager = AMP::Discretization::DOFManager::shared_ptr( 
      new AMP::Discretization::multiDOFManager( subchannelMesh->getComm(), subchannelChildrenDOFManagers ) );
  }

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  (new AMP::LinearAlgebra::Variable("flow"));
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable (new AMP::LinearAlgebra::Variable("flow"));

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr SolVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr ResVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get nonlinear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelFourEqNonlinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = subchannelMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
  subchannelOpParams->d_dofMap = subchannelDOFManager;
  subchannelOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
  subchannelOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
  subchannelOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");
  // create nonlinear operator
  boost::shared_ptr<AMP::Operator::SubchannelFourEqNonlinearOperator> subchannelOperator (new AMP::Operator::SubchannelFourEqNonlinearOperator(subchannelOpParams));
  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);

  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // check number of lateral gaps
  std::map<std::vector<double>,AMP::Mesh::MeshElement> lateralFaceMap = subchannelOperator->getLateralFaces(subchannelOpParams->d_Mesh);
  size_t Ngaps = lateralFaceMap.size();
  if (Ngaps == 108) {// for 3x3 subchannel array with 9 axial intervals, there are 12x9=108 gaps
     ut->passes(exeName+": number of lateral gaps");
  } else {
     std::cout<<"Incorrent number of lateral gaps. Found: "<<Ngaps<<". Expected: 108."<<std::endl;
     ut->failure(exeName+": number of lateral gaps");
  }

  // number of subchannels
  const size_t numSubchannels = 3*3; // 3x3 subchannel array
  // number of axial intervals
  const size_t numAxialIntervals = 9;
  // compute height of subchannels
  std::vector<double> box = subchannelOpParams->d_Mesh->getBoundingBox();
  const double height = box[5] - box[4];
  // height of each axial interval
  const double dz = height/numAxialIntervals;
  // get all of the unique x,y,z points in subchannel mesh
  subchannelOperator->fillSubchannelGrid(subchannelOpParams->d_Mesh);

  // put all cells in an array by subchannel
  AMP::Mesh::MeshElement d_elem[numSubchannels][numAxialIntervals]; // array of array of elements for each subchannel
  AMP::Mesh::MeshIterator cell = subchannelOpParams->d_Mesh->getIterator(AMP::Mesh::Volume, 0); // iterator for cells of mesh
  for( ; cell != cell.end(); ++cell) { // loop over all cells
    std::vector<double> center = cell->centroid();
    // get the index of the subchannel
    int isub = subchannelOperator->getSubchannelIndex( center[0], center[1] );
    bool found_axial_index = false;
    // get the axial index of the cell
    for (unsigned int j = 0; j < numAxialIntervals; ++j)
    {
       // center of interval j
       double center_interval_j = (j+0.5)*dz;
       // check if center of cell is equal to center of axial interval j
       if (AMP::Utilities::approx_equal(center[2],center_interval_j,1.0e-12))
       {
          d_elem[isub][j] = *cell;
          found_axial_index = true;
          break;
       }
    }
    // ensure that a valid axial index was found
    if (!found_axial_index) AMP_ERROR("A cell was not in center of any axial interval");
  }// end for cell

  // get scale factors for axial mass flow rate, enthalpy, and pressure
  const double m_scale = AMP::Operator::Subchannel::scaleAxialMassFlowRate;
  const double h_scale = AMP::Operator::Subchannel::scaleEnthalpy;
  const double p_scale = AMP::Operator::Subchannel::scalePressure;

  // set axial face quantities of input vector for the nonlinear residual
  // loop over subchannels
  for(size_t isub =0; isub < numSubchannels; ++isub){
      // loop over axial intervals
      for (unsigned int j = 0; j < numAxialIntervals; ++j) {
         // get axial faces
         AMP::Mesh::MeshElement plusFace;  // upper axial face for current cell
         AMP::Mesh::MeshElement minusFace; // lower axial face for current cell
         subchannelOperator->getAxialFaces(d_elem[isub][j],plusFace,minusFace);

         if (j == 0) // for first axial interval only, set the quantities for the lower face
         {
            // get dofs on minus face
            std::vector<size_t> minusDofs;
            subchannelDOFManager->getDOFs(minusFace.globalID(),minusDofs);
            // set values of minus face
            SolVec->setValueByGlobalID(minusDofs[0], m_scale*0.35);
            SolVec->setValueByGlobalID(minusDofs[1], h_scale*1000.0e3);
            SolVec->setValueByGlobalID(minusDofs[2], p_scale*15.5e6);
         }

         // get dofs on plus face
         std::vector<size_t> plusDofs;
         subchannelDOFManager->getDOFs(plusFace.globalID(),plusDofs);
         // set values of plus face
         SolVec->setValueByGlobalID(plusDofs[0], m_scale*0.35);
         SolVec->setValueByGlobalID(plusDofs[1], h_scale*1000.0e3);
         SolVec->setValueByGlobalID(plusDofs[2], p_scale*15.5e6);
      }
  }

  // get scale factors for lateral mass flow rate
  const double w_scale = AMP::Operator::Subchannel::scaleLateralMassFlowRate;

  // set lateral face quantities (lateral mass flow rates) of input vector for nonlinear residual
  // loop over all faces in mesh
  AMP::Mesh::MeshIterator face = subchannelOpParams->d_Mesh->getIterator(AMP::Mesh::Face, 0);
  for (; face != face.end(); face++) { // loop over all faces in mesh
     std::vector<double> faceCentroid = face->centroid();
     // try to find face in lateral face map
     std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
     if (lateralFaceIterator != lateralFaceMap.end()) { // if face in lateral face map,
        // get lateral face
        AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
        // get crossflow from solution vector
        std::vector<size_t> gapDofs;
        subchannelDOFManager->getDOFs(lateralFace.globalID(),gapDofs);
        SolVec->setValueByGlobalID(gapDofs[0], w_scale*0.001);
     }
  }

  // apply the operator
  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

  // initialize success boolean for known residual comparison test
  bool passedKnownTest = true;
  // subchannel 1: check inlet and outlet residuals
  size_t isub = 0;
  // extract subchannel cells from d_elem[isub]
/*
  boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > subchannelElements( new std::vector<AMP::Mesh::MeshElement>() );
  subchannelElements->reserve(numSubchannels);
  for(size_t ielem=0; ielem < d_elem[isub].size(); ++ielem){
    subchannelElements->push_back(d_elem[isub][ielem]);
  }
  AMP::Mesh::MeshIterator     localSubchannelCell = AMP::Mesh::MultiVectorIterator( subchannelElements ); // iterator over elements of current subchannel
*/

  double m_inlet_residual = 0.;
  double h_inlet_residual = 0.;
  double p_inlet_residual = 0.;
  double m_outlet_residual = 0.;
  double h_outlet_residual = 0.;
  double p_outlet_residual = 0.;
  // loop over axial intervals
  for (size_t j = 0; j < numAxialIntervals; ++j) {
     // get axial faces
     AMP::Mesh::MeshElement plusFace;  // upper axial face for current cell
     AMP::Mesh::MeshElement minusFace; // lower axial face for current cell
     subchannelOperator->getAxialFaces(d_elem[isub][j],plusFace,minusFace);

     // get unknowns of first face
     if (j == 0)
     {
        // get dofs on minus face
        std::vector<size_t> minusDofs;
        subchannelDOFManager->getDOFs(minusFace.globalID(),minusDofs);
        // get values of minus face
        m_inlet_residual = ResVec->getValueByGlobalID(minusDofs[0])/m_scale;
        h_inlet_residual = ResVec->getValueByGlobalID(minusDofs[1])/h_scale;
        p_inlet_residual = ResVec->getValueByGlobalID(minusDofs[2])/p_scale;
     }
    
     // get unknowns on last face
     if (j == numAxialIntervals-1)
     {
        // get dofs on plus face
        std::vector<size_t> plusDofs;
        subchannelDOFManager->getDOFs(plusFace.globalID(),plusDofs);
        // get values of plus face
        m_outlet_residual = ResVec->getValueByGlobalID(plusDofs[0])/m_scale;
        h_outlet_residual = ResVec->getValueByGlobalID(plusDofs[1])/h_scale;
        p_outlet_residual = ResVec->getValueByGlobalID(plusDofs[2])/p_scale;
     }
  }

  // compare residual with MATLAB residual values
  if (!AMP::Utilities::approx_equal(m_inlet_residual,0.038,1.0e-6)) passedKnownTest = false;
  if (!AMP::Utilities::approx_equal(h_inlet_residual,-2.630791e5,1.0e-6)) passedKnownTest = false;
  if (!AMP::Utilities::approx_equal(p_inlet_residual,3.782389e-1,1.0e-6)) passedKnownTest = false;
  if (!AMP::Utilities::approx_equal(m_outlet_residual,0.002,1.0e-6)) passedKnownTest = false;
  if (!AMP::Utilities::approx_equal(h_outlet_residual,1.185235e3,1.0e-6)) passedKnownTest = false;
  if (!AMP::Utilities::approx_equal(p_outlet_residual,0.0,1.0e-6)) passedKnownTest = false;

  if (passedKnownTest) ut->passes(exeName+": known value test");
  else ut->failure(exeName+": known residual test");

}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);

    AMP::UnitTest ut;

    // create list of input files
    const int NUMFILES=1;
    std::string files[NUMFILES] = {
        "testSubchannelFourEqNonlinearOperator"
    };

    // execute unit test for each input file
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


