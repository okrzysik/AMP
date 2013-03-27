#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <iomanip>
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

void compare_face_value(std::string variable, unsigned int i, unsigned int j, double value_array[][10], double known_value, bool &passed_known_test)
{
   if (!AMP::Utilities::approx_equal(value_array[i][j],known_value,1.0e-6)) {
      std::cout << "Residual value for " << variable << "[" << i << "][" << j << "]" << " does not match known MATLAB value" << std::endl;
      passed_known_test = false;
   }
}

void compare_gap_value(std::string variable, unsigned int i, unsigned int j, double value_array[][9], double known_value, bool &passed_known_test)
{
   if (!AMP::Utilities::approx_equal(value_array[i][j],known_value,1.0e-6)) {
      std::cout << "Residual value for " << variable << "[" << i << "][" << j << "]" << " does not match known MATLAB value" << std::endl;
      passed_known_test = false;
   }
}

unsigned int getMATLABGapIndex(AMP::Mesh::MeshElement gapFace)
{
   double pitch = 0.0126; // pitch for test problem [m]
   double x1 = 0.5*pitch;
   double x2 = 1.0*pitch;
   double x3 = 1.5*pitch;
   double x4 = 2.0*pitch;
   double x5 = 2.5*pitch;

   // get gap face centroid
   std::vector<double> centroid = gapFace.centroid();
   // gap MATLAB index
   unsigned int k;
   // look at location of gap to determine gap MATLAB index
   if      ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 1;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 2;
   else if ((AMP::Utilities::approx_equal(centroid[0],x5,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 3;
   else if ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 4;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 5;
   else if ((AMP::Utilities::approx_equal(centroid[0],x5,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 6;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x5,1.0e-12))) k = 7;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) k = 8;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) k = 9;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x5,1.0e-12))) k = 10;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) k = 11;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) k = 12;
   else AMP_ERROR("Gap face did not match any known location");

   return k;
}

unsigned int getMATLABAxialIndex(AMP::Mesh::MeshElement gapFace)
{
   double height = 3.66;
   unsigned int Nz = 9;
   double dz = height/Nz;

   // get gap face centroid
   std::vector<double> centroid = gapFace.centroid();
   // axial interval MATLAB index
   unsigned int j;
   // boolean for if the axial index has been found
   bool foundIndex = false;
   // loop over axial intervals
   for (unsigned int i = 0; i < Nz; ++i)
   {
      if (AMP::Utilities::approx_equal(centroid[2],(i+0.5)*dz,1.0e-12))
      {
         j = i + 1;
         foundIndex = true;
         break;
      }
   }
   if (!foundIndex) AMP_ERROR("Axial index was not found for gap face");
  
   return j;
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
  // number of gaps
  const size_t numGaps = 12;
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
      size_t ii = isub+1; // corresponding MATLAB index for this subchannel
      // loop over axial intervals
      for (unsigned int j = 0; j < numAxialIntervals; ++j) {
         // get axial faces
         AMP::Mesh::MeshElement plusFace;  // upper axial face for current cell
         AMP::Mesh::MeshElement minusFace; // lower axial face for current cell
         subchannelOperator->getAxialFaces(d_elem[isub][j],plusFace,minusFace);

         if (j == 0) // for first axial interval only, set the quantities for the lower face
         {
            size_t jj = j+1; // corresponding MATLAB index for this axial face
            // get dofs on minus face
            std::vector<size_t> minusDofs;
            subchannelDOFManager->getDOFs(minusFace.globalID(),minusDofs);
            // set values of minus face
            SolVec->setValueByGlobalID(minusDofs[0], m_scale*0.35     *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
            SolVec->setValueByGlobalID(minusDofs[1], h_scale*1000.0e3 *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
            SolVec->setValueByGlobalID(minusDofs[2], p_scale*15.5e6   *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
         }

         size_t jj = j+2; // corresponding MATLAB index for this axial face
         // get dofs on plus face
         std::vector<size_t> plusDofs;
         subchannelDOFManager->getDOFs(plusFace.globalID(),plusDofs);
         // set values of plus face
         SolVec->setValueByGlobalID(plusDofs[0], m_scale*0.35     *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
         SolVec->setValueByGlobalID(plusDofs[1], h_scale*1000.0e3 *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
         SolVec->setValueByGlobalID(plusDofs[2], p_scale*15.5e6   *(1.0 + 1.0/100.0*cos(ii)*cos(17.3*jj)));
      }
  }

  // get scale factors for lateral mass flow rate
  const double w_scale = AMP::Operator::Subchannel::scaleLateralMassFlowRate;

  // set lateral face quantities (lateral mass flow rates) of input vector for nonlinear residual
  // array of gap faces
  AMP::Mesh::MeshElement gapFaces[numGaps][numAxialIntervals]; // gap faces
  // loop over all faces in mesh
  AMP::Mesh::MeshIterator face = subchannelOpParams->d_Mesh->getIterator(AMP::Mesh::Face, 0);
  for (; face != face.end(); face++) { // loop over all faces in mesh
     std::vector<double> faceCentroid = face->centroid();
     // try to find face in lateral face map
     std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = lateralFaceMap.find(faceCentroid);
     if (lateralFaceIterator != lateralFaceMap.end()) { // if face in lateral face map,
        // get lateral face
        AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
        // get MATLAB index for gap
        unsigned int k = getMATLABGapIndex(lateralFace);
        // get MATLAB axial index
        unsigned int j = getMATLABAxialIndex(lateralFace);
        // put gap face in array
        gapFaces[k-1][j-1] = lateralFace;
        // get crossflow from solution vector
        std::vector<size_t> gapDofs;
        subchannelDOFManager->getDOFs(lateralFace.globalID(),gapDofs);
        // set test value for crossflow
        SolVec->setValueByGlobalID(gapDofs[0], w_scale*0.001*(1.0 + 1.0/100.0*cos(k)*cos(17.3*j)));
     }
  }

  // apply the operator
  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

  // initialize success boolean for known residual comparison test
  bool passed_known_test = true;

  std::cout << std::setprecision(6) << std::scientific;

  double m_res[numSubchannels][numAxialIntervals+1];
  double h_res[numSubchannels][numAxialIntervals+1];
  double p_res[numSubchannels][numAxialIntervals+1];
  // loop over subchannels
  for (size_t isub = 0; isub < numSubchannels; ++isub) {
     size_t ii = isub+1; // corresponding MATLAB index for this subchannel
     std::cout << "\nSubchannel " << ii << ":" << std::endl;
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
           m_res[ii-1][0] = ResVec->getValueByGlobalID(minusDofs[0])/m_scale;
           h_res[ii-1][0] = ResVec->getValueByGlobalID(minusDofs[1])/h_scale;
           p_res[ii-1][0] = ResVec->getValueByGlobalID(minusDofs[2])/p_scale;
           std::cout << "Face 0:\t" << m_res[ii-1][0] << "\t" << h_res[ii-1][0] << "\t" << p_res[ii-1][0] << std::endl;
        }
       
        // get dofs on plus face
        std::vector<size_t> plusDofs;
        subchannelDOFManager->getDOFs(plusFace.globalID(),plusDofs);
        // get values of plus face
        m_res[ii-1][j+1] = ResVec->getValueByGlobalID(plusDofs[0])/m_scale;
        h_res[ii-1][j+1] = ResVec->getValueByGlobalID(plusDofs[1])/h_scale;
        p_res[ii-1][j+1] = ResVec->getValueByGlobalID(plusDofs[2])/p_scale;
        std::cout << "Face " << j+1 << ":\t" << m_res[ii-1][j+1] << "\t" << h_res[ii-1][j+1] << "\t" << p_res[ii-1][j+1] << std::endl;
     }
  }

  double w_res[numGaps][numAxialIntervals];
  // loop over gaps
  for (size_t k = 0; k < numGaps; ++k) {
     std::cout << "\nGap " << k+1 << ":" << std::endl;
     for (size_t j = 0; j < numAxialIntervals; ++j) {
        std::vector<size_t> gapDofs;
        subchannelDOFManager->getDOFs(gapFaces[k][j].globalID(),gapDofs);
        w_res[k][j] = ResVec->getValueByGlobalID(gapDofs[0])/w_scale;
        std::cout << "Interval " << j+1 << ":\t" << w_res[k][j] << std::endl;
     }
  }

  // compare residual with MATLAB residual values
  compare_face_value("Axial mass flow rate",0,0,m_res,0.0380401638196053,passed_known_test);
  compare_face_value("Axial mass flow rate",0,1,m_res,7.07590468003155e-05,passed_known_test);
  compare_face_value("Axial mass flow rate",0,2,m_res,0.00375600264827943,passed_known_test);
  compare_face_value("Axial mass flow rate",0,3,m_res,0.00400383175999823,passed_known_test);
  compare_face_value("Axial mass flow rate",0,4,m_res,0.000329115348446149,passed_known_test);
  compare_face_value("Axial mass flow rate",0,5,m_res,-7.48069571874154e-05,passed_known_test);
  compare_face_value("Axial mass flow rate",0,6,m_res,0.00358275179763167,passed_known_test);
  compare_face_value("Axial mass flow rate",0,7,m_res,0.00414203848082026,passed_known_test);
  compare_face_value("Axial mass flow rate",0,8,m_res,0.000508236891157803,passed_known_test);
  compare_face_value("Axial mass flow rate",0,9,m_res,-0.000205405021838726,passed_known_test);
  compare_face_value("Axial mass flow rate",1,0,m_res,0.0379690653800828,passed_known_test);
  compare_face_value("Axial mass flow rate",1,1,m_res,0.0024857078056411,passed_known_test);
  compare_face_value("Axial mass flow rate",1,2,m_res,-0.000342378775625912,passed_known_test);
  compare_face_value("Axial mass flow rate",1,3,m_res,-0.00054272885626487,passed_known_test);
  compare_face_value("Axial mass flow rate",1,4,m_res,0.00227684733114156,passed_known_test);
  compare_face_value("Axial mass flow rate",1,5,m_res,0.00259696628733119,passed_known_test);
  compare_face_value("Axial mass flow rate",1,6,m_res,-0.000209012009868793,passed_known_test);
  compare_face_value("Axial mass flow rate",1,7,m_res,-0.000648322235645622,passed_known_test);
  compare_face_value("Axial mass flow rate",1,8,m_res,0.00213899521032629,passed_known_test);
  compare_face_value("Axial mass flow rate",1,9,m_res,0.00269670403721126,passed_known_test);
  compare_face_value("Axial mass flow rate",2,0,m_res,0.0379264080874497,passed_known_test);
  compare_face_value("Axial mass flow rate",2,1,m_res,0.00353540758780944,passed_known_test);
  compare_face_value("Axial mass flow rate",2,2,m_res,-0.00323969692872978,passed_known_test);
  compare_face_value("Axial mass flow rate",2,3,m_res,-0.00367302219102955,passed_known_test);
  compare_face_value("Axial mass flow rate",2,4,m_res,0.00308367570126416,passed_known_test);
  compare_face_value("Axial mass flow rate",2,5,m_res,0.00380400938440724,passed_known_test);
  compare_face_value("Axial mass flow rate",2,6,m_res,-0.0029220904500517,passed_known_test);
  compare_face_value("Axial mass flow rate",2,7,m_res,-0.00392813282147675,passed_known_test);
  compare_face_value("Axial mass flow rate",2,8,m_res,0.00275523273109986,passed_known_test);
  compare_face_value("Axial mass flow rate",2,9,m_res,0.00404516854037586,passed_known_test);
  compare_face_value("Axial mass flow rate",3,0,m_res,0.0379514108598289,passed_known_test);
  compare_face_value("Axial mass flow rate",3,1,m_res,0.00333399337998071,passed_known_test);
  compare_face_value("Axial mass flow rate",3,2,m_res,-0.00112662659789425,passed_known_test);
  compare_face_value("Axial mass flow rate",3,3,m_res,-0.00142432740872004,passed_known_test);
  compare_face_value("Axial mass flow rate",3,4,m_res,0.00302364694832642,passed_known_test);
  compare_face_value("Axial mass flow rate",3,5,m_res,0.00351028710741343,passed_known_test);
  compare_face_value("Axial mass flow rate",3,6,m_res,-0.000917015935656216,passed_known_test);
  compare_face_value("Axial mass flow rate",3,7,m_res,-0.00159171737485982,passed_known_test);
  compare_face_value("Axial mass flow rate",3,8,m_res,0.00280692595932617,passed_known_test);
  compare_face_value("Axial mass flow rate",3,9,m_res,0.00366847128252814,passed_known_test);
  compare_face_value("Axial mass flow rate",4,0,m_res,0.0380210862636011,passed_known_test);
  compare_face_value("Axial mass flow rate",4,1,m_res,-0.00101282774497257,passed_known_test);
  compare_face_value("Axial mass flow rate",4,2,m_res,0.000920211529605241,passed_known_test);
  compare_face_value("Axial mass flow rate",4,3,m_res,0.00105191613728551,passed_known_test);
  compare_face_value("Axial mass flow rate",4,4,m_res,-0.000875528639496079,passed_known_test);
  compare_face_value("Axial mass flow rate",4,5,m_res,-0.00108910650693511,passed_known_test);
  compare_face_value("Axial mass flow rate",4,6,m_res,0.000829265991133191,passed_known_test);
  compare_face_value("Axial mass flow rate",4,7,m_res,0.00112433174955635,passed_known_test);
  compare_face_value("Axial mass flow rate",4,8,m_res,-0.000781507058428509,passed_known_test);
  compare_face_value("Axial mass flow rate",4,9,m_res,-0.00115752830655695,passed_known_test);
  compare_face_value("Axial mass flow rate",5,0,m_res,0.0380713750538627,passed_known_test);
  compare_face_value("Axial mass flow rate",5,1,m_res,-0.00442852596361711,passed_known_test);
  compare_face_value("Axial mass flow rate",5,2,m_res,0.00212412796804168,passed_known_test);
  compare_face_value("Axial mass flow rate",5,3,m_res,0.00256123147286351,passed_known_test);
  compare_face_value("Axial mass flow rate",5,4,m_res,-0.00397285534381356,passed_known_test);
  compare_face_value("Axial mass flow rate",5,5,m_res,-0.0046875112816139,passed_known_test);
  compare_face_value("Axial mass flow rate",5,6,m_res,0.00181621865431677,passed_known_test);
  compare_face_value("Axial mass flow rate",5,7,m_res,0.00280713753716213,passed_known_test);
  compare_face_value("Axial mass flow rate",5,8,m_res,-0.00365450052663179,passed_known_test);
  compare_face_value("Axial mass flow rate",5,9,m_res,-0.00491989439212611,passed_known_test);
  compare_face_value("Axial mass flow rate",6,0,m_res,0.0380560419487658,passed_known_test);
  compare_face_value("Axial mass flow rate",6,1,m_res,-0.00269237399766216,passed_known_test);
  compare_face_value("Axial mass flow rate",6,2,m_res,0.0024708251781098,passed_known_test);
  compare_face_value("Axial mass flow rate",6,3,m_res,0.00279732876163915,passed_known_test);
  compare_face_value("Axial mass flow rate",6,4,m_res,-0.00235200132017567,passed_known_test);
  compare_face_value("Axial mass flow rate",6,5,m_res,-0.00289723617136585,passed_known_test);
  compare_face_value("Axial mass flow rate",6,6,m_res,0.00222893363362082,passed_known_test);
  compare_face_value("Axial mass flow rate",6,7,m_res,0.00299191595911848,passed_known_test);
  compare_face_value("Axial mass flow rate",6,8,m_res,-0.00210184417536605,passed_known_test);
  compare_face_value("Axial mass flow rate",6,9,m_res,-0.00308119728962148,passed_known_test);
  compare_face_value("Axial mass flow rate",7,0,m_res,0.0379891841344244,passed_known_test);
  compare_face_value("Axial mass flow rate",7,1,m_res,-0.000480080957017993,passed_known_test);
  compare_face_value("Axial mass flow rate",7,2,m_res,-0.00149106255146655,passed_known_test);
  compare_face_value("Axial mass flow rate",7,3,m_res,-0.00154077820964579,passed_known_test);
  compare_face_value("Axial mass flow rate",7,4,m_res,-0.000531908417851741,passed_known_test);
  compare_face_value("Axial mass flow rate",7,5,m_res,-0.000439338375712318,passed_known_test);
  compare_face_value("Axial mass flow rate",7,6,m_res,-0.00144427601277116,passed_known_test);
  compare_face_value("Axial mass flow rate",7,7,m_res,-0.00157953341031048,passed_known_test);
  compare_face_value("Axial mass flow rate",7,8,m_res,-0.00058034118509241,passed_known_test);
  compare_face_value("Axial mass flow rate",7,9,m_res,-0.000402640483553133,passed_known_test);
  compare_face_value("Axial mass flow rate",8,0,m_res,0.0379322703770131,passed_known_test);
  compare_face_value("Axial mass flow rate",8,1,m_res,0.00125342539339376,passed_known_test);
  compare_face_value("Axial mass flow rate",8,2,m_res,-0.00496498829548058,passed_known_test);
  compare_face_value("Axial mass flow rate",8,3,m_res,-0.00537937102868342,passed_known_test);
  compare_face_value("Axial mass flow rate",8,4,m_res,0.0008214406694291,passed_known_test);
  compare_face_value("Axial mass flow rate",8,5,m_res,0.00149921910298596,passed_known_test);
  compare_face_value("Axial mass flow rate",8,6,m_res,-0.00467280218285727,passed_known_test);
  compare_face_value("Axial mass flow rate",8,7,m_res,-0.00561275336868163,passed_known_test);
  compare_face_value("Axial mass flow rate",8,8,m_res,0.000519341031305047,passed_known_test);
  compare_face_value("Axial mass flow rate",8,9,m_res,0.00171976897045735,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,0,w_res,0.0010001147537703,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,1,w_res,7523.0278487063,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,2,w_res,479.484740864891,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,3,w_res,-7502.6604740481,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,4,w_res,-798.17986790002,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,5,w_res,7468.75570163437,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,6,w_res,1115.4348133507,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,7,w_res,-7421.37468378016,passed_known_test);
  compare_gap_value("Lateral mass flow rate",0,8,w_res,-1430.67711377483,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,0,w_res,0.000999911615371665,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,1,w_res,-4409.48173524556,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,2,w_res,-281.040997967558,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,3,w_res,4397.54378947444,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,4,w_res,467.838181963105,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,5,w_res,-4377.67113274599,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,6,w_res,-653.791210802806,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,7,w_res,4349.89964644907,passed_known_test);
  compare_gap_value("Lateral mass flow rate",1,8,w_res,838.564585161023,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,0,w_res,0.000999789737392714,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,1,w_res,-12287.9341446666,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,2,w_res,-783.178936253159,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,3,w_res,12254.6665757808,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,4,w_res,1303.72796797571,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,5,w_res,-12199.2873138265,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,6,w_res,-1821.9246077612,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,7,w_res,12121.8963048463,passed_known_test);
  compare_gap_value("Lateral mass flow rate",2,8,w_res,2336.83387482741,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,0,w_res,0.000999861173885226,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,1,w_res,-8868.91656728573,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,2,w_res,-565.265769652848,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,3,w_res,8844.90543050937,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,4,w_res,940.976275393172,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,5,w_res,-8804.93499556926,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,6,w_res,-1314.98891986421,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,7,w_res,8749.07740645474,passed_known_test);
  compare_gap_value("Lateral mass flow rate",3,8,w_res,1686.62887955874,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,0,w_res,0.00100006024646743,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,1,w_res,2704.14200364906,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,2,w_res,172.350141626911,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,3,w_res,-2696.82097450005,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,4,w_res,-286.904662339548,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,5,w_res,2684.63395428302,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,6,w_res,400.941519479304,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,7,w_res,-2667.60290828019,passed_known_test);
  compare_gap_value("Lateral mass flow rate",4,8,w_res,-514.254926389997,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,0,w_res,0.00100020392872532,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,1,w_res,11791.0248896967,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,2,w_res,751.508130665646,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,3,w_res,-11759.1026101011,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,4,w_res,-1251.00677351256,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,5,w_res,11705.9628298768,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,6,w_res,1748.24817797759,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,7,w_res,-11631.7014089308,passed_known_test);
  compare_gap_value("Lateral mass flow rate",5,8,w_res,-2242.33512150492,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,0,w_res,0.00100016011985362,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,1,w_res,6026.56566159104,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,2,w_res,384.106838592368,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,3,w_res,-6010.24971279669,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,4,w_res,-639.407913958829,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,5,w_res,5983.08919461331,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,6,w_res,893.555287417077,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,7,w_res,-5945.13309004638,passed_known_test);
  compare_gap_value("Lateral mass flow rate",6,8,w_res,-1146.09036387998,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,0,w_res,0.000999969097526927,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,1,w_res,-5905.94393896087,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,2,w_res,-376.418913235937,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,3,w_res,5889.95456663249,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,4,w_res,626.610149720149,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,5,w_res,-5863.33765523439,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,6,w_res,-875.670751059293,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,7,w_res,5826.14125488916,passed_known_test);
  compare_gap_value("Lateral mass flow rate",7,8,w_res,1123.15135015794,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,0,w_res,0.000999806486791466,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,1,w_res,5667.1147417398,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,2,w_res,361.196953407278,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,3,w_res,-5651.77195034251,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,4,w_res,-601.270747810733,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,5,w_res,5626.231407954,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,6,w_res,840.259653250444,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,7,w_res,-5590.53917503865,passed_known_test);
  compare_gap_value("Lateral mass flow rate",8,8,w_res,-1077.73242526632,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,0,w_res,0.0009998217912075,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,1,w_res,3615.78929456195,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,2,w_res,230.454430227205,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,3,w_res,-3606.00013390738,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,4,w_res,-383.628714322006,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,5,w_res,3589.70450640959,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,6,w_res,536.110811219483,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,7,w_res,-3566.93179113409,passed_known_test);
  compare_gap_value("Lateral mass flow rate",9,8,w_res,-687.625565374865,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,0,w_res,0.00100000093996551,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,1,w_res,-4262.66310964657,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,2,w_res,-271.683361594061,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,3,w_res,4251.12265376985,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,4,w_res,452.260924955529,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,5,w_res,-4231.91168198128,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,6,w_res,-632.022439574428,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,7,w_res,4205.06488100564,passed_known_test);
  compare_gap_value("Lateral mass flow rate",10,8,w_res,810.643578522305,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,0,w_res,0.00100017922452356,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,1,w_res,4824.21970461761,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,2,w_res,307.47456187076,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,3,w_res,-4811.15891245038,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,4,w_res,-511.841118052397,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,5,w_res,4789.41712745467,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,6,w_res,715.28414764942,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,7,w_res,-4759.03355624499,passed_known_test);
  compare_gap_value("Lateral mass flow rate",11,8,w_res,-917.436542865777,passed_known_test);
  compare_face_value("Specific Enthalpy",0,0,h_res,-262962.972855096,passed_known_test);
  compare_face_value("Specific Enthalpy",0,1,h_res,-3310.27497242607,passed_known_test);
  compare_face_value("Specific Enthalpy",0,2,h_res,762.111541951088,passed_known_test);
  compare_face_value("Specific Enthalpy",0,3,h_res,-510.447982107353,passed_known_test);
  compare_face_value("Specific Enthalpy",0,4,h_res,-8822.47793126787,passed_known_test);
  compare_face_value("Specific Enthalpy",0,5,h_res,-10567.0810023223,passed_known_test);
  compare_face_value("Specific Enthalpy",0,6,h_res,-3309.23291756892,passed_known_test);
  compare_face_value("Specific Enthalpy",0,7,h_res,-276.666575083377,passed_known_test);
  compare_face_value("Specific Enthalpy",0,8,h_res,-4740.38763875071,passed_known_test);
  compare_face_value("Specific Enthalpy",0,9,h_res,-3779.01055660154,passed_known_test);
  compare_face_value("Specific Enthalpy",1,0,h_res,-263168.515839132,passed_known_test);
  compare_face_value("Specific Enthalpy",1,1,h_res,2319.92253390948,passed_known_test);
  compare_face_value("Specific Enthalpy",1,2,h_res,-6202.61678797473,passed_known_test);
  compare_face_value("Specific Enthalpy",1,3,h_res,-9352.8977824601,passed_known_test);
  compare_face_value("Specific Enthalpy",1,4,h_res,-5639.62422912767,passed_known_test);
  compare_face_value("Specific Enthalpy",1,5,h_res,-5341.56023838203,passed_known_test);
  compare_face_value("Specific Enthalpy",1,6,h_res,-10121.3070907309,passed_known_test);
  compare_face_value("Specific Enthalpy",1,7,h_res,-9543.4814485272,passed_known_test);
  compare_face_value("Specific Enthalpy",1,8,h_res,-1728.60489068715,passed_known_test);
  compare_face_value("Specific Enthalpy",1,9,h_res,2699.91563067949,passed_known_test);
  compare_face_value("Specific Enthalpy",2,0,h_res,-263291.836531722,passed_known_test);
  compare_face_value("Specific Enthalpy",2,1,h_res,5574.88742622351,passed_known_test);
  compare_face_value("Specific Enthalpy",2,2,h_res,-10289.487958886,passed_known_test);
  compare_face_value("Specific Enthalpy",2,3,h_res,-13968.9185471552,passed_known_test);
  compare_face_value("Specific Enthalpy",2,4,h_res,-2669.46596768912,passed_known_test);
  compare_face_value("Specific Enthalpy",2,5,h_res,-1160.72707033023,passed_known_test);
  compare_face_value("Specific Enthalpy",2,6,h_res,-13499.9503695949,passed_known_test);
  compare_face_value("Specific Enthalpy",2,7,h_res,-14430.3782960175,passed_known_test);
  compare_face_value("Specific Enthalpy",2,8,h_res,524.20568961014,passed_known_test);
  compare_face_value("Specific Enthalpy",2,9,h_res,6491.47681203384,passed_known_test);
  compare_face_value("Specific Enthalpy",3,0,h_res,-263219.554417629,passed_known_test);
  compare_face_value("Specific Enthalpy",3,1,h_res,4012.92664487938,passed_known_test);
  compare_face_value("Specific Enthalpy",3,2,h_res,-7194.85437812648,passed_known_test);
  compare_face_value("Specific Enthalpy",3,3,h_res,-11070.0823288619,passed_known_test);
  compare_face_value("Specific Enthalpy",3,4,h_res,-4737.16034625536,passed_known_test);
  compare_face_value("Specific Enthalpy",3,5,h_res,-3568.34936012506,passed_known_test);
  compare_face_value("Specific Enthalpy",3,6,h_res,-10964.9432749866,passed_known_test);
  compare_face_value("Specific Enthalpy",3,7,h_res,-11335.2355947579,passed_known_test);
  compare_face_value("Specific Enthalpy",3,8,h_res,-977.111628063604,passed_known_test);
  compare_face_value("Specific Enthalpy",3,9,h_res,4540.36307251752,passed_known_test);
  compare_face_value("Specific Enthalpy",4,0,h_res,-263018.125366054,passed_known_test);
  compare_face_value("Specific Enthalpy",4,1,h_res,-3992.56105099566,passed_known_test);
  compare_face_value("Specific Enthalpy",4,2,h_res,-4187.01136354822,passed_known_test);
  compare_face_value("Specific Enthalpy",4,3,h_res,-6635.06753263771,passed_known_test);
  compare_face_value("Specific Enthalpy",4,4,h_res,-12120.5350580285,passed_known_test);
  compare_face_value("Specific Enthalpy",4,5,h_res,-13524.0499477232,passed_known_test);
  compare_face_value("Specific Enthalpy",4,6,h_res,-9371.86754430794,passed_known_test);
  compare_face_value("Specific Enthalpy",4,7,h_res,-6518.37932080872,passed_known_test);
  compare_face_value("Specific Enthalpy",4,8,h_res,-6930.58513504395,passed_known_test);
  compare_face_value("Specific Enthalpy",4,9,h_res,-4225.99506463542,passed_known_test);
  compare_face_value("Specific Enthalpy",5,0,h_res,-262872.742334926,passed_known_test);
  compare_face_value("Specific Enthalpy",5,1,h_res,-9474.90026811425,passed_known_test);
  compare_face_value("Specific Enthalpy",5,2,h_res,-645.248834960898,passed_known_test);
  compare_face_value("Specific Enthalpy",5,3,h_res,-1289.60778952278,passed_known_test);
  compare_face_value("Specific Enthalpy",5,4,h_res,-14982.6630097227,passed_known_test);
  compare_face_value("Specific Enthalpy",5,5,h_res,-17862.9318018612,passed_known_test);
  compare_face_value("Specific Enthalpy",5,6,h_res,-5486.8142384355,passed_known_test);
  compare_face_value("Specific Enthalpy",5,7,h_res,-887.573805167401,passed_known_test);
  compare_face_value("Specific Enthalpy",5,8,h_res,-10121.3690754197,passed_known_test);
  compare_face_value("Specific Enthalpy",5,9,h_res,-10283.4275264019,passed_known_test);
  compare_face_value("Specific Enthalpy",6,0,h_res,-262917.06977172,passed_known_test);
  compare_face_value("Specific Enthalpy",6,1,h_res,-6855.67106696425,passed_known_test);
  compare_face_value("Specific Enthalpy",6,2,h_res,40.7384773168261,passed_known_test);
  compare_face_value("Specific Enthalpy",6,3,h_res,-1049.56439343089,passed_known_test);
  compare_face_value("Specific Enthalpy",6,4,h_res,-12300.9901124822,passed_known_test);
  compare_face_value("Specific Enthalpy",6,5,h_res,-14369.1827403579,passed_known_test);
  compare_face_value("Specific Enthalpy",6,6,h_res,-4242.98678594535,passed_known_test);
  compare_face_value("Specific Enthalpy",6,7,h_res,-708.105037486514,passed_known_test);
  compare_face_value("Specific Enthalpy",6,8,h_res,-8001.11126945905,passed_known_test);
  compare_face_value("Specific Enthalpy",6,9,h_res,-7541.19901847183,passed_known_test);
  compare_face_value("Specific Enthalpy",7,0,h_res,-263110.353252648,passed_known_test);
  compare_face_value("Specific Enthalpy",7,1,h_res,-1620.07113410543,passed_known_test);
  compare_face_value("Specific Enthalpy",7,2,h_res,-6597.62001661329,passed_known_test);
  compare_face_value("Specific Enthalpy",7,3,h_res,-9381.88936306201,passed_known_test);
  compare_face_value("Specific Enthalpy",7,4,h_res,-9191.08048473349,passed_known_test);
  compare_face_value("Specific Enthalpy",7,5,h_res,-9438.69034441606,passed_known_test);
  compare_face_value("Specific Enthalpy",7,6,h_res,-10701.6675370741,passed_known_test);
  compare_face_value("Specific Enthalpy",7,7,h_res,-9446.03252752563,passed_known_test);
  compare_face_value("Specific Enthalpy",7,8,h_res,-5089.80350539571,passed_known_test);
  compare_face_value("Specific Enthalpy",7,9,h_res,-1491.81079242963,passed_known_test);
  compare_face_value("Specific Enthalpy",8,0,h_res,-263274.888863299,passed_known_test);
  compare_face_value("Specific Enthalpy",8,1,h_res,3030.38830830804,passed_known_test);
  compare_face_value("Specific Enthalpy",8,2,h_res,-11650.9588302808,passed_known_test);
  compare_face_value("Specific Enthalpy",8,3,h_res,-15275.9125590784,passed_known_test);
  compare_face_value("Specific Enthalpy",8,4,h_res,-5066.13351359758,passed_known_test);
  compare_face_value("Specific Enthalpy",8,5,h_res,-3632.9679914306,passed_known_test);
  compare_face_value("Specific Enthalpy",8,6,h_res,-14849.3233611213,passed_known_test);
  compare_face_value("Specific Enthalpy",8,7,h_res,-15692.6724177729,passed_known_test);
  compare_face_value("Specific Enthalpy",8,8,h_res,-1883.08726075716,passed_known_test);
  compare_face_value("Specific Enthalpy",8,9,h_res,3858.49223320375,passed_known_test);
  compare_face_value("Pressure",0,0,p_res,-4.37247372415232,passed_known_test);
  compare_face_value("Pressure",0,1,p_res,10.6191939888063,passed_known_test);
  compare_face_value("Pressure",0,2,p_res,11.5843585873826,passed_known_test);
  compare_face_value("Pressure",0,3,p_res,-2.80084106462981,passed_known_test);
  compare_face_value("Pressure",0,4,p_res,-4.3768306130786,passed_known_test);
  compare_face_value("Pressure",0,5,p_res,9.94120985653401,passed_known_test);
  compare_face_value("Pressure",0,6,p_res,12.1256010129246,passed_known_test);
  compare_face_value("Pressure",0,7,p_res,-2.0998472144618,passed_known_test);
  compare_face_value("Pressure",0,8,p_res,-5.45475985548001,passed_known_test);
  compare_face_value("Pressure",0,9,p_res,-81864.8092918955,passed_known_test);
  compare_face_value("Pressure",1,0,p_res,8.9813893184275,passed_known_test);
  compare_face_value("Pressure",1,1,p_res,-1.55899565477114,passed_known_test);
  compare_face_value("Pressure",1,2,p_res,-2.301478110143,passed_known_test);
  compare_face_value("Pressure",1,3,p_res,8.77820592847964,passed_known_test);
  compare_face_value("Pressure",1,4,p_res,9.99146908665285,passed_known_test);
  compare_face_value("Pressure",1,5,p_res,-1.03681542626484,passed_known_test);
  compare_face_value("Pressure",1,6,p_res,-2.71840220086999,passed_known_test);
  compare_face_value("Pressure",1,7,p_res,8.23832651523694,passed_known_test);
  compare_face_value("Pressure",1,8,p_res,9.81503725306052,passed_known_test);
  compare_face_value("Pressure",1,9,p_res,63053.185302658,passed_known_test);
  compare_face_value("Pressure",2,0,p_res,16.9914887879986,passed_known_test);
  compare_face_value("Pressure",2,1,p_res,-8.86824512599643,passed_known_test);
  compare_face_value("Pressure",2,2,p_res,-10.6339253762245,passed_known_test);
  compare_face_value("Pressure",2,3,p_res,15.7241985504114,passed_known_test);
  compare_face_value("Pressure",2,4,p_res,18.6102532355529,passed_known_test);
  compare_face_value("Pressure",2,5,p_res,-7.62599119290171,passed_known_test);
  compare_face_value("Pressure",2,6,p_res,-11.6258180037208,passed_known_test);
  compare_face_value("Pressure",2,7,p_res,14.4398758029329,passed_known_test);
  compare_face_value("Pressure",2,8,p_res,18.9747253574555,passed_known_test);
  compare_face_value("Pressure",2,9,p_res,150000.37211461,passed_known_test);
  compare_face_value("Pressure",3,0,p_res,12.2994693136461,passed_known_test);
  compare_face_value("Pressure",3,1,p_res,-4.5803355496898,passed_known_test);
  compare_face_value("Pressure",3,2,p_res,-5.74901875049736,passed_known_test);
  compare_face_value("Pressure",3,3,p_res,11.653105902382,passed_known_test);
  compare_face_value("Pressure",3,4,p_res,13.5613146242567,passed_known_test);
  compare_face_value("Pressure",3,5,p_res,-3.76006514189337,passed_known_test);
  compare_face_value("Pressure",3,6,p_res,-6.40374755756331,passed_known_test);
  compare_face_value("Pressure",3,7,p_res,10.805052235449,passed_known_test);
  compare_face_value("Pressure",3,8,p_res,13.6085918290409,passed_known_test);
  compare_face_value("Pressure",3,9,p_res,99037.9085665476,passed_known_test);
  compare_face_value("Pressure",4,0,p_res,-0.797664973468099,passed_known_test);
  compare_face_value("Pressure",4,1,p_res,7.34307212098663,passed_known_test);
  compare_face_value("Pressure",4,2,p_res,7.85021117811943,passed_known_test);
  compare_face_value("Pressure",4,3,p_res,0.298056632375483,passed_known_test);
  compare_face_value("Pressure",4,4,p_res,-0.529804531305929,passed_known_test);
  compare_face_value("Pressure",4,5,p_res,6.98711112988511,passed_known_test);
  compare_face_value("Pressure",4,6,p_res,8.13434403411397,passed_known_test);
  compare_face_value("Pressure",4,7,p_res,0.666091370348251,passed_known_test);
  compare_face_value("Pressure",4,8,p_res,-1.36581605845869,passed_known_test);
  compare_face_value("Pressure",4,9,p_res,-42979.5513808839,passed_known_test);
  compare_face_value("Pressure",5,0,p_res,-10.2512460320998,passed_known_test);
  compare_face_value("Pressure",5,1,p_res,15.9481002837834,passed_known_test);
  compare_face_value("Pressure",5,2,p_res,17.6643785897697,passed_known_test);
  compare_face_value("Pressure",5,3,p_res,-7.89940058356422,passed_known_test);
  compare_face_value("Pressure",5,4,p_res,-10.7008670246342,passed_known_test);
  compare_face_value("Pressure",5,5,p_res,14.7432414055283,passed_known_test);
  compare_face_value("Pressure",5,6,p_res,18.6261537079572,passed_known_test);
  compare_face_value("Pressure",5,7,p_res,-6.65362717425516,passed_known_test);
  compare_face_value("Pressure",5,8,p_res,-12.1744930108092,passed_known_test);
  compare_face_value("Pressure",5,9,p_res,-145481.809999086,passed_known_test);
  compare_face_value("Pressure",6,0,p_res,-7.36490302198328,passed_known_test);
  compare_face_value("Pressure",6,1,p_res,13.3289753236269,passed_known_test);
  compare_face_value("Pressure",6,2,p_res,14.6746793352881,passed_known_test);
  compare_face_value("Pressure",6,3,p_res,-5.39795822277163,passed_known_test);
  compare_face_value("Pressure",6,4,p_res,-7.59586909072042,passed_known_test);
  compare_face_value("Pressure",6,5,p_res,12.3829935107838,passed_known_test);
  compare_face_value("Pressure",6,6,p_res,15.4299530396493,passed_known_test);
  compare_face_value("Pressure",6,7,p_res,-4.41985674561837,passed_known_test);
  compare_face_value("Pressure",6,8,p_res,-8.87519327381994,passed_known_test);
  compare_face_value("Pressure",6,9,p_res,-114228.76342787,passed_known_test);
  compare_face_value("Pressure",7,0,p_res,5.19195092440458,passed_known_test);
  compare_face_value("Pressure",7,1,p_res,1.87648979681475,passed_known_test);
  compare_face_value("Pressure",7,2,p_res,1.61667274719758,passed_known_test);
  compare_face_value("Pressure",7,3,p_res,5.49055381954823,passed_known_test);
  compare_face_value("Pressure",7,4,p_res,5.91495546933588,passed_known_test);
  compare_face_value("Pressure",7,5,p_res,2.05907078996018,passed_known_test);
  compare_face_value("Pressure",7,6,p_res,1.47091128343418,passed_known_test);
  compare_face_value("Pressure",7,7,p_res,5.30178171739982,passed_known_test);
  compare_face_value("Pressure",7,8,p_res,5.48340999498541,passed_known_test);
  compare_face_value("Pressure",7,9,p_res,22045.6814459972,passed_known_test);
  compare_face_value("Pressure",8,0,p_res,15.8806787496428,passed_known_test);
  compare_face_value("Pressure",8,1,p_res,-7.87354706143448,passed_known_test);
  compare_face_value("Pressure",8,2,p_res,-9.49913873690652,passed_known_test);
  compare_face_value("Pressure",8,3,p_res,14.7591065544117,passed_known_test);
  compare_face_value("Pressure",8,4,p_res,17.4157629192827,passed_known_test);
  compare_face_value("Pressure",8,5,p_res,-6.73024002259814,passed_known_test);
  compare_face_value("Pressure",8,6,p_res,-10.4119822273427,passed_known_test);
  compare_face_value("Pressure",8,7,p_res,13.5770786117365,passed_known_test);
  compare_face_value("Pressure",8,8,p_res,17.705869688255,passed_known_test);
  compare_face_value("Pressure",8,9,p_res,138051.428467283,passed_known_test);

  if (passed_known_test) ut->passes(exeName+": known value test");
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


