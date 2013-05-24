#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

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
#include "operators/subchannel/SubchannelFourEqLinearOperator.h"
#include "operators/subchannel/SubchannelConstants.h"
#include "operators/OperatorBuilder.h"

#include "ampmesh/StructuredMeshHelper.h"
#include "discretization/simpleDOF_Manager.h"
#include "discretization/MultiDOF_Manager.h"

// number of subchannels
const size_t numSubchannels = 3*3; // 3x3 subchannel array
// number of axial intervals
const size_t numAxialIntervals = 3;
// number of axial faces
const size_t numAxialFaces = numAxialIntervals + 1;
// number of gaps in MATLAB
const size_t numGaps_MATLAB = 12;
// total number of DoFs in MATLAB
const unsigned int num_dofs_MATLAB = 3*numSubchannels*numAxialFaces + numGaps_MATLAB*numAxialIntervals;
// number of gaps in AMP; AMP automatically puts DOFs on exterior gaps, as opposed to MATLAB
const size_t numGaps_AMP = 24;
// total number of DoFs in AMP
const unsigned int num_dofs_AMP = 3*numSubchannels*numAxialFaces + numGaps_AMP*numAxialIntervals;

// function to get the MATLAB subchannel index
size_t getMATLABSubchannelIndex(AMP::Mesh::MeshElement face)
{
   double pitch = 0.0126; // pitch for test problem [m]
   double x1 = 0.5*pitch;
   double x2 = 1.5*pitch;
   double x3 = 2.5*pitch;

   // get face centroid
   std::vector<double> centroid = face.centroid();
   // gap MATLAB index
   size_t i;
   // look at location of subchannel to determine subchannel MATLAB index
   if      ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) i = 0;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) i = 1;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) i = 2;
   else if ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) i = 3;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) i = 4;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) i = 5;
   else if ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) i = 6;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) i = 7;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) i = 8;
   else AMP_ERROR("Face did not match any known location");

   return i;
}

// function to get the MATLAB gap index
size_t getMATLABGapIndex(AMP::Mesh::MeshElement gapFace)
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
   size_t k = 0;
   // look at location of gap to determine gap MATLAB index
   if      ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 0;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 1;
   else if ((AMP::Utilities::approx_equal(centroid[0],x5,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x4,1.0e-12))) k = 2;
   else if ((AMP::Utilities::approx_equal(centroid[0],x1,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 3;
   else if ((AMP::Utilities::approx_equal(centroid[0],x3,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 4;
   else if ((AMP::Utilities::approx_equal(centroid[0],x5,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x2,1.0e-12))) k = 5;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x5,1.0e-12))) k = 6;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) k = 7;
   else if ((AMP::Utilities::approx_equal(centroid[0],x2,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) k = 8;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x5,1.0e-12))) k = 9;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x3,1.0e-12))) k = 10;
   else if ((AMP::Utilities::approx_equal(centroid[0],x4,1.0e-12))&&(AMP::Utilities::approx_equal(centroid[1],x1,1.0e-12))) k = 11;
   else AMP_ERROR("Gap face did not match any known location");

   return k;
}

// function to determine the axial level
size_t getMATLABAxialIndex(AMP::Mesh::MeshElement face, bool is_axial_face_quantity)
{
   double height = 3.66;
   double dz = height / numAxialIntervals;
   std::vector<double> centroid = face.centroid();
   size_t j = 0;
   // boolean for if the axial index has been found
   bool foundIndex = false;
   // loop over axial intervals
   if (is_axial_face_quantity) {
      for (size_t k = 0; k < numAxialFaces; ++k)
      {
         if (AMP::Utilities::approx_equal(centroid[2],k*dz,1.0e-12))
         {
            j = k + 1;
            foundIndex = true;
            break;
         }
      }
   } else {
      for (size_t k = 0; k < numAxialIntervals; ++k)
      {
         if (AMP::Utilities::approx_equal(centroid[2],(k+0.5)*dz,1.0e-12))
         {
            j = k + 1;
            foundIndex = true;
            break;
         }
      }
   }
   if (!foundIndex) AMP_ERROR("Axial index was not found");

   return j;
}

// function used to get all lateral gaps
std::map<std::vector<double>,AMP::Mesh::MeshElement> getLateralFaces(
   AMP::Mesh::Mesh::shared_ptr mesh,
   bool getExterior
)
{
   // map of lateral gaps to their centroids
   std::map<std::vector<double>,AMP::Mesh::MeshElement> lateralFaceMap;
   // get iterator over all faces of mesh
   AMP::Mesh::MeshIterator face = mesh->getIterator(AMP::Mesh::Face,0);
   // loop over faces
   for (; face != face.end(); face++) {
      // check that face is vertical
      // ---------------------------
      // get centroid of current face
      std::vector<double> faceCentroid = face->centroid();
      // get vertices of current face
      std::vector<AMP::Mesh::MeshElement> vertices = face->getElements(AMP::Mesh::Vertex);

      bool perpindicular_to_x = true; // is the current face perpindicular to x-axis?
      bool perpindicular_to_y = true; // is the current face perpindicular to y-axis?
      // loop over vertices of current face
      for (size_t j=0; j<vertices.size(); ++j) {
        // get coordinates of current vertex
        std::vector<double> vertexCoord = vertices[j].coord();
        // if any vertex does not have the same x-coordinate as the face centroid,
        if ( !AMP::Utilities::approx_equal(vertexCoord[0],faceCentroid[0], 1.0e-6) )
          // then the face is not perpindicular to x-axis
          perpindicular_to_x = false;
        // if any vertex does not have the same y-coordinate as the face centroid,
        if ( !AMP::Utilities::approx_equal(vertexCoord[1],faceCentroid[1], 1.0e-6) )
          // then the face is not perpindicular to y-axis
          perpindicular_to_y = false;
      }
      // check that face is in the interior of the mesh; it must have two adjacent cells
      // -------------------------------------------------------------------------------
      // if the face is vertical
      if (perpindicular_to_x || perpindicular_to_y) {
         // if the face has more than 1 adjacent cell
         if ((mesh->getElementParents(*face,AMP::Mesh::Volume)).size() > 1) {
            // insert face into map with centroid
            lateralFaceMap.insert(std::pair<std::vector<double>,AMP::Mesh::MeshElement>(faceCentroid,*face));
         }
      }
   }// end loop over faces
   return lateralFaceMap;
}

// find the MATLAB index for a variable on a given face
size_t AMP_to_MATLAB(AMP::Mesh::MeshElement face, size_t variable_id) {
   // determine if variable is axial face quantity or lateral face quantity
   bool is_axial_face_quantity;
   switch (variable_id) {
      case 0: {
         is_axial_face_quantity = true;
         break;
      } case 1: {
         is_axial_face_quantity = true;
         break;
      } case 2: {
         is_axial_face_quantity = true;
         break;
      } case 3: {
         is_axial_face_quantity = false;
         break;
      } default: {
         AMP_ERROR("Invalid variable id");
         break;
      }
   }

   // determine axial level
   size_t j = getMATLABAxialIndex(face,is_axial_face_quantity);

   // determine subchannel or gap
   size_t i;
   if (is_axial_face_quantity) i = getMATLABSubchannelIndex(face);
   else                        i = getMATLABGapIndex(face);

   // compute index
   size_t index;
   switch (variable_id) {
      case 0: {
         index = i*numAxialFaces + j;
         break;
      } case 1: {
         index = i*numAxialFaces + j + numGaps_MATLAB*numAxialIntervals + numSubchannels*numAxialFaces;
         break;
      } case 2: {
         index = i*numAxialFaces + j + numGaps_MATLAB*numAxialIntervals + 2*numSubchannels*numAxialFaces;
         break;
      } case 3: {
         index = i*numAxialIntervals + j + numSubchannels*numAxialFaces;
         break;
      } default: {
         AMP_ERROR("Invalid variable id");
         break;
      }
   }

   return index;
}

// function to create map of global IDs to elements and variables
void createGlobalIDMaps(
   AMP::Discretization::DOFManager::shared_ptr dof_manager,
   boost::shared_ptr<AMP::Mesh::Mesh>          mesh,
   std::map<size_t,AMP::Mesh::MeshElement>     &elements_by_globalID,
   std::map<size_t,size_t>                     &variables_by_globalID
)
{
   // loop over axial faces  
   AMP::Mesh::Mesh::shared_ptr xyMesh =
      mesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(mesh, 0));
   AMP::Mesh::MeshIterator axial_face = xyMesh->getIterator(AMP::Mesh::Face,0);
   for (; axial_face != axial_face.end(); axial_face++) {
      std::vector<size_t> dofs;
      dof_manager->getDOFs(axial_face->globalID(),dofs);
      elements_by_globalID.insert(std::make_pair(dofs[0],*axial_face));
      elements_by_globalID.insert(std::make_pair(dofs[1],*axial_face));
      elements_by_globalID.insert(std::make_pair(dofs[2],*axial_face));
      variables_by_globalID.insert(std::make_pair(dofs[0],0));
      variables_by_globalID.insert(std::make_pair(dofs[1],1));
      variables_by_globalID.insert(std::make_pair(dofs[2],2));
   }
  
   // loop over lateral faces
   std::map<std::vector<double>,AMP::Mesh::MeshElement> lateral_face_map =
      getLateralFaces(mesh,true);
   AMP::Mesh::MeshIterator face = mesh->getIterator(AMP::Mesh::Face,0);
   for (; face != face.end(); face++) {
      std::vector<double> centroid = face->centroid();
      std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateral_face_iterator =
         lateral_face_map.find(centroid);
      if (lateral_face_iterator != lateral_face_map.end()) {
         AMP::Mesh::MeshElement lateral_face = lateral_face_iterator->second;
         std::vector<size_t> dofs;
         dof_manager->getDOFs(lateral_face.globalID(),dofs);
         elements_by_globalID.insert(std::make_pair(dofs[0],lateral_face));
         variables_by_globalID.insert(std::make_pair(dofs[0],3));
      }
   }
}

// function to check that Jacobian matches known values
bool JacobianIsCorrect(
   boost::shared_ptr<AMP::LinearAlgebra::Matrix> J_test_AMP,
   double                                        J_reference[num_dofs_MATLAB][num_dofs_MATLAB],
   AMP::Discretization::DOFManager::shared_ptr   dof_manager,
   AMP::Mesh::Mesh::shared_ptr                   mesh,
   std::map<std::vector<double>,AMP::Mesh::MeshElement> lateral_face_map,
   std::map<size_t,AMP::Mesh::MeshElement>       elements_by_globalID,
   std::map<size_t,size_t>                       variables_by_globalID
)
{
   double J_test_MATLAB[num_dofs_MATLAB][num_dofs_MATLAB] = {0.0};
   // loop over axial faces  
   AMP::Mesh::Mesh::shared_ptr xyMesh =
      mesh->Subset( AMP::Mesh::StructuredMeshHelper::getXYFaceIterator(mesh, 0));
   AMP::Mesh::MeshIterator axial_face = xyMesh->getIterator(AMP::Mesh::Face,0);
   for (; axial_face != axial_face.end(); axial_face++) {
      std::vector<size_t> dofs;
      dof_manager->getDOFs(axial_face->globalID(),dofs);
      size_t i_m_AMP = dofs[0];
      size_t i_h_AMP = dofs[1];
      size_t i_p_AMP = dofs[2];
      size_t i_m_MATLAB = AMP_to_MATLAB(*axial_face,0);
      size_t i_h_MATLAB = AMP_to_MATLAB(*axial_face,1);
      size_t i_p_MATLAB = AMP_to_MATLAB(*axial_face,2);
      std::vector<double> val_m, val_h, val_p; // nonzero values in Jacobian row
      std::vector<unsigned int> ind_m, ind_h, ind_p; // indices of nonzeros in Jacobian row
      J_test_AMP->getRowByGlobalID(i_m_AMP, ind_m, val_m);
      J_test_AMP->getRowByGlobalID(i_h_AMP, ind_h, val_h);
      J_test_AMP->getRowByGlobalID(i_p_AMP, ind_p, val_p);
      // loop over all DOFs
      for (size_t j_AMP = 0; j_AMP < num_dofs_AMP; j_AMP++) {
         std::map<size_t,AMP::Mesh::MeshElement>::iterator face_iterator = elements_by_globalID.find(j_AMP);
         if (face_iterator != elements_by_globalID.end()) {
            AMP::Mesh::MeshElement face = face_iterator->second;
            std::map<size_t,size_t>::iterator variable_id_iterator = variables_by_globalID.find(j_AMP);
            if (variable_id_iterator == variables_by_globalID.end()) AMP_ERROR("index exists only for face");
            size_t variable_id = variable_id_iterator->second;
            size_t j_MATLAB = AMP_to_MATLAB(face,variable_id);
            std::vector<unsigned int>::iterator nonzero_ind_m = std::find(ind_m.begin(),ind_m.end(),j_AMP);
            std::vector<unsigned int>::iterator nonzero_ind_h = std::find(ind_h.begin(),ind_h.end(),j_AMP);
            std::vector<unsigned int>::iterator nonzero_ind_p = std::find(ind_p.begin(),ind_p.end(),j_AMP);
            if (nonzero_ind_m != ind_m.end()) J_test_MATLAB[i_m_MATLAB][j_MATLAB] = val_m[*nonzero_ind_m];
            if (nonzero_ind_h != ind_h.end()) J_test_MATLAB[i_h_MATLAB][j_MATLAB] = val_h[*nonzero_ind_h];
            if (nonzero_ind_p != ind_p.end()) J_test_MATLAB[i_p_MATLAB][j_MATLAB] = val_p[*nonzero_ind_p];
         }
      }
   }
   
   // loop over lateral faces
   AMP::Mesh::MeshIterator lateral_face = mesh->getIterator(AMP::Mesh::Face,0);
   for (; lateral_face != lateral_face.end(); lateral_face++) {
      std::vector<double> centroid = lateral_face->centroid();
      std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateral_face_iterator =
         lateral_face_map.find(centroid);
      if (lateral_face_iterator != lateral_face_map.end()) {
         std::vector<size_t> dofs;
         dof_manager->getDOFs(lateral_face->globalID(),dofs);
         size_t i_w_AMP = dofs[0];
         size_t i_w_MATLAB = AMP_to_MATLAB(*lateral_face,3);
         std::vector<double> val_w; // nonzero values in Jacobian row
         std::vector<unsigned int> ind_w; // indices of nonzeros in Jacobian row
         J_test_AMP->getRowByGlobalID(i_w_AMP, ind_w, val_w);
         // loop over all DOFs
         for (size_t j_AMP = 0; j_AMP < num_dofs_AMP; j_AMP++) {
            std::map<size_t,AMP::Mesh::MeshElement>::iterator face_iterator = elements_by_globalID.find(j_AMP);
            if (face_iterator != elements_by_globalID.end()) {
               AMP::Mesh::MeshElement face = face_iterator->second;
               std::map<size_t,size_t>::iterator variable_id_iterator = variables_by_globalID.find(j_AMP);
               if (variable_id_iterator == variables_by_globalID.end()) AMP_ERROR("index exists only for face");
               size_t variable_id = variable_id_iterator->second;
               size_t j_MATLAB = AMP_to_MATLAB(face,variable_id);
               std::vector<unsigned int>::iterator nonzero_ind_w = std::find(ind_w.begin(),ind_w.end(),j_AMP);
               if (nonzero_ind_w != ind_w.end()) J_test_MATLAB[i_w_MATLAB][j_MATLAB] = val_w[*nonzero_ind_w];
            }
         }
      }
   }

   // compare each matrix entry
   bool passed = true;
   for (size_t i = 0; i < num_dofs_MATLAB; i++) {
      for (size_t j = 0; j < num_dofs_MATLAB; j++) {
         if (!AMP::Utilities::approx_equal(J_reference[i][j],J_test_MATLAB[i][j],1.0e-6)) {
            passed = false;
            std::cout << "Values are inconsistent for entry (" << i << "," << j << ") (MATLAB indices): " <<
               J_reference[i][j] << " (MATLAB) vs. " << J_test_MATLAB[i][j] << " (AMP)" << std::endl;
         }
      }
   }

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

  // check number of cells
  size_t expected_number_of_cells = numSubchannels * numAxialIntervals;
  size_t numCells = subchannelMesh->numGlobalElements(AMP::Mesh::Volume);
  if (numCells == expected_number_of_cells) ut->passes(exeName+": number of cells");
  else {                                    ut->failure(exeName+": number of cells");
     std::cout << "number of cells in mesh: " << numCells << std::endl;
     std::cout << "expected number of cells: " << expected_number_of_cells << std::endl;
  }

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

  // check number of DOFs
  if (subchannelDOFManager->numGlobalDOF() == num_dofs_AMP) ut->passes(exeName+": number of DOFs");
  else {                                                ut->failure(exeName+": number of DOFs");
     std::cout << "number of DOFs in DOF Manager: " << subchannelDOFManager->numGlobalDOF() << std::endl;
     std::cout << "expected number of DOFs: " << num_dofs_AMP << std::endl;
  }

  // get input and output variables
  AMP::LinearAlgebra::Variable::shared_ptr inputVariable  (new AMP::LinearAlgebra::Variable("flow"));
  AMP::LinearAlgebra::Variable::shared_ptr outputVariable (new AMP::LinearAlgebra::Variable("flow"));

  // create solution, rhs, and residual vectors
  AMP::LinearAlgebra::Vector::shared_ptr FrozenVec = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr SolVec    = AMP::LinearAlgebra::createVector( subchannelDOFManager, inputVariable , true );
  AMP::LinearAlgebra::Vector::shared_ptr RhsVec    = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );
  AMP::LinearAlgebra::Vector::shared_ptr ResVec    = AMP::LinearAlgebra::createVector( subchannelDOFManager, outputVariable, true );

  // create subchannel physics model
  boost::shared_ptr<AMP::Database> subchannelPhysics_db = input_db->getDatabase("SubchannelPhysicsModel");
  boost::shared_ptr<AMP::Operator::ElementPhysicsModelParameters> params( new AMP::Operator::ElementPhysicsModelParameters(subchannelPhysics_db));
  boost::shared_ptr<AMP::Operator::SubchannelPhysicsModel>  subchannelPhysicsModel (new AMP::Operator::SubchannelPhysicsModel(params));

  // get linear operator database
  boost::shared_ptr<AMP::Database> subchannelOperator_db = input_db->getDatabase("SubchannelFourEqLinearOperator");
  // set operator parameters
  boost::shared_ptr<AMP::Operator::SubchannelOperatorParameters> subchannelOpParams(new AMP::Operator::SubchannelOperatorParameters( subchannelOperator_db ));
  subchannelOpParams->d_Mesh = subchannelMesh ;
  subchannelOpParams->d_subchannelPhysicsModel = subchannelPhysicsModel;
  subchannelOpParams->d_dofMap = subchannelDOFManager;
  subchannelOpParams->clad_x = input_db->getDatabase("CladProperties")->getDoubleArray("x");
  subchannelOpParams->clad_y = input_db->getDatabase("CladProperties")->getDoubleArray("y");
  subchannelOpParams->clad_d = input_db->getDatabase("CladProperties")->getDoubleArray("d");
  // create linear operator
  boost::shared_ptr<AMP::Operator::SubchannelFourEqLinearOperator> subchannelOperator (new AMP::Operator::SubchannelFourEqLinearOperator(subchannelOpParams));
  // reset the nonlinear operator
  subchannelOperator->reset(subchannelOpParams);

  // report successful creation
  ut->passes(exeName+": creation");
  std::cout.flush();

  // check number of lateral gaps
  std::map<std::vector<double>,AMP::Mesh::MeshElement> interiorLateralFaceMap = getLateralFaces(subchannelMesh,false);
  size_t Ngaps = interiorLateralFaceMap.size();
  if (Ngaps == 36) {// for 3x3 subchannel array with 3 axial intervals, there are 12x3=36 gaps
     ut->passes(exeName+": number of lateral gaps");
  } else {
     std::cout<<"Incorrent number of lateral gaps. Found: "<<Ngaps<<". Expected: 36."<<std::endl;
     ut->failure(exeName+": number of lateral gaps");
  }

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

  // set axial face quantities of input vector
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

  // set lateral face quantities (lateral mass flow rates) of input vector
  // array of gap faces
  AMP::Mesh::MeshElement gapFaces[numGaps_MATLAB][numAxialIntervals]; // gap faces
  // loop over all faces in mesh
  AMP::Mesh::MeshIterator face = subchannelOpParams->d_Mesh->getIterator(AMP::Mesh::Face, 0);
  for (; face != face.end(); face++) { // loop over all faces in mesh
     std::vector<double> faceCentroid = face->centroid();
     // try to find face in lateral face map
     std::map<std::vector<double>,AMP::Mesh::MeshElement>::iterator lateralFaceIterator = interiorLateralFaceMap.find(faceCentroid);
     if (lateralFaceIterator != interiorLateralFaceMap.end()) { // if face in lateral face map,
        // get lateral face
        AMP::Mesh::MeshElement lateralFace = lateralFaceIterator->second;
        // get MATLAB index for gap
        unsigned int k = getMATLABGapIndex(lateralFace);
        // get MATLAB axial index
        unsigned int j = getMATLABAxialIndex(lateralFace,false);
        // put gap face in array
        gapFaces[k][j] = lateralFace;
        // get crossflow from solution vector
        std::vector<size_t> gapDofs;
        subchannelDOFManager->getDOFs(lateralFace.globalID(),gapDofs);
        // set test value for crossflow
        SolVec->setValueByGlobalID(gapDofs[0], w_scale*0.001*(1.0 + 1.0/100.0*cos(k)*cos(17.3*j)));
     }
  }

  // apply the operator
  subchannelOperator->setFrozenVector(FrozenVec);
  subchannelOperator->reset(subchannelOpParams);
//  subchannelOperator->apply(RhsVec, SolVec, ResVec, 1.0, 0.0);

  // get the AMP Jacobian matrix to be tested against the MATLAB Jacobian matrix
  boost::shared_ptr<AMP::LinearAlgebra::Matrix> testJacobian = subchannelOperator->getMatrix();

  // get the MATLAB Jacobian matrix
  double knownJacobian[num_dofs_MATLAB][num_dofs_MATLAB] = {
    {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-0.00677668901258425,0.00674966917997032,0,0,0,0,0,0,0,0,0,0,-0.00680134179247133,0.00675123052041697,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75450738336373,4.75720888252606,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-0.00674966917997032,0.00679583334404178,0,0,0,0,0,0,0,0,0,0,-0.00675123052041697,0.00677129943182299,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75113294598402,4.7536054258804,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-0.00679502285863625,0.00681578079453035,0,0,0,0,0,0,0,0,0,0,-0.00678057423704236,0.00681485684336233,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75504157219458,4.75730493425882,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-0.00681578079453035,0.00678027402447358,0,0,0,0,0,0,0,0,0,0,-0.00681485684336233,0.00679466062218694,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75117047528547,4.75296033111564,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,-0.00680615306651423,0.00685545819708779,0,0,0,0,0,0,0,0,0,0,-0.0067659218851098,0.00685286883884361,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75207026445087,4.75713741578303,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,-0.00685545819708779,0.00677105865990493,0,0,0,0,0,0,0,0,0,0,-0.00685286883884361,0.00681119707306303,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75096878063051,4.75633680932167,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00679961730689246,0.00683220093186155,0,0,0,0,0,0,0,0,0,0,-0.00677060673094952,0.006830338373976,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74793016128515,4.75685852673517,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00683220093186155,0.00677644930146954,0,0,0,0,0,0,0,0,0,0,-0.006830338373976,0.00680544001532756,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75071169563154,4.76023829062968,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00678158213030341,0.00676740617221729,0,0,0,0,0,0,0,0,0,0,-0.00679042650059103,0.00676796877440433,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74664515776453,4.75672569526946,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00676740617221729,0.0067916335529895,0,0,0,0,0,0,0,0,0,0,-0.00676796877440433,0.00678279735427316,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75063646484676,4.76129322816635,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00676872486067267,0.00672065503545349,0,0,0,0,0,0,0,0,0,0,-0.00680732255528395,0.00672309165630821,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74961025253147,4.75686201176696,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00672065503545349,0.00680274407584255,0,0,0,0,0,0,0,0,0,0,-0.00672309165630821,0.00676425457951943,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75081312459266,4.75874514045713,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0.366193650793651,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-0.00677699640914768,0.00673518671850398,0,0,-0.0067967116015305,0.00673643465839852,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74712413718295,4.75672497663301,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-0.00673518671850398,0.0067949086811636,0,0,-0.00673643465839852,0.00677522361754851,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75066428209911,4.76090435637508,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00680035124824462,0.00679773991114397,0,0,-0.00678096403351543,0.00679650559033574,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.76242481837535,4.75775930419679,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00679773991114397,0.00677425751608255,0,0,-0.00679650559033574,0.00679348992193275,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75163913917039,4.74565951727907,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00677023641473485,0.00684779571028473,0,0,-0.0067887034679692,0.00684898861682152,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.73215050067433,4.75586677521402,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00684779571028473,0.00680655538876296,0,0,-0.00684898861682152,0.00678795878231737,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74970902130792,4.77587467194847,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-0.00679441244743652,0.0068444611640022,0,0,-0.00680637127556524,0.00684522244914879,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.79457858836221,4.75985586012582,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-0.0068444611640022,0.00678209953915717,0,0,-0.00684522244914879,0.00677040834669981,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.7536836042023,4.71375709187597,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00678117996261892,0.00678633985029537,0,0,-0.00676735116574352,0.00678545097362945,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.71253324828324,4.75455171703698,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00678633985029537,0.00679284140848072,0,0,-0.00678545097362945,0.00680688675891539,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74845565171368,4.79553782624389,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00679123436366546,0.00672998171707532,0,0,-0.00680715442010169,0.00673098033097224,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.78372216201276,4.7590600383179,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.00672998171707532,0.00678037736113117,0,0,-0.00673098033097224,0.00676475590104976,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75299393772574,4.7245020804524,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.366193650793651,0,0,0,0.366193650793651,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1000094.19706307,994622.408115175,0,0,9.14366044264196,9.14366044264196,0,0,0,0,0,0,11.413523529005,11.413523529005,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000114.7537703,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000114.7537703,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.20800981525427,0.348110647992782,0,0,-0.0699966637489306,0,0,0,0,0,0,0,-0.0700334099427809,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-995569.394593763,998688.402559005,0,0,-429.944638472618,-429.944638472618,0,0,0,0,0,0,-536.55066884904,-536.55066884904,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,994601.851407948,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,994601.851407948,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.206093443437389,0.349879581010687,0,0,-0.0699989268532521,0,0,0,0,0,0,0,-0.0700312080717666,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-999717.581831606,1005321.89790061,0,0,-27.4118739980478,-27.4118739980478,0,0,0,0,0,0,-34.2199460998972,-34.2199460998972,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,999655.94574482,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,999655.94574482,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.207847040069733,0.351884236895589,0,0,-0.0700343676693115,0,0,0,0,0,0,0,-0.0699989973965461,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-9.14347408824019,-9.14347408824019,0,0,-999921.964878841,1004147.36448996,0,0,5.48397791432918,5.48397791432918,0,0,0,0,0,0,-6.68965780523463,-6.68965780523463,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,999911.615371665,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000114.7537703,0,0,0,0,0,0,0,0,999911.615371665,0,0,0,0,0,0,0,0,-0.0709968238687842,0,0,0,-0.137702903320502,0.351455199898998,0,0,-0.0702333947721765,0,0,0,0,0,0,0,-0.0700363701318951,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,430.354641697736,430.354641697736,0,0,-1003671.07843432,1000751.6299712,0,0,-257.725050945294,-257.725050945294,0,0,0,0,0,0,314.78117189742,314.78117189742,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1004157.71399714,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-994601.851407948,0,0,0,0,0,0,0,0,1004157.71399714,0,0,0,0,0,0,0,0,-0.0709913946322184,0,0,0,-0.139202699538987,0.350092748042935,0,0,-0.0702071651132457,0,0,0,0,0,0,0,-0.0700338675341108,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,27.4135418549702,27.4135418549702,0,0,-1000233.97716082,995884.559581242,0,0,-16.4518000525577,-16.4518000525577,0,0,0,0,0,0,20.0586665873575,20.0586665873575,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000264.99440839,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999655.94574482,0,0,0,0,0,0,0,0,1000264.99440839,0,0,0,0,0,0,0,0,-0.0710338875986626,0,0,0,-0.138255021256138,0.348548739816786,0,0,-0.0698064182258784,0,0,0,0,0,0,0,-0.0699961415923722,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-5.48391083553225,-5.48391083553225,0,0,-999813.864463332,1009866.86644102,0,0,0,0,0,0,0,0,0,0,-18.6423179641248,-18.6423179641248,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,999789.737392714,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999911.615371665,0,0,0,0,0,0,0,0,0,0,0,0,-0.071233216563384,0,0,0,-0.208672731442553,0.353461847729074,0,0,0,0,0,0,0,0,0,0,-0.070020492135328,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,257.871135015198,257.871135015198,0,0,-1008757.15657443,1001764.24537217,0,0,0,0,0,0,0,0,0,0,877.814667261428,877.814667261428,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1009890.99351164,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1004157.71399714,0,0,0,0,0,0,0,0,0,0,0,0,-0.0712155482586517,0,0,0,-0.212225085865979,0.350220642952239,0,0,0,0,0,0,0,0,0,0,-0.0700197057563376,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,16.4524001963693,16.4524001963693,0,0,-1000558.06489177,990208.128279203,0,0,0,0,0,0,0,0,0,0,55.8986767779937,55.8986767779937,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000630.40843497,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000264.99440839,0,0,0,0,0,0,0,0,0,0,0,0,-0.0708069525307052,0,0,0,-0.209402145907097,0.346547524657601,0,0,0,0,0,0,0,0,0,0,-0.0700114484108286,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-11.413233138775,-11.413233138775,0,0,0,0,0,0,0,0,0,0,-999895.003154797,1006496.70993662,0,0,-8.95942964008051,-8.95942964008051,0,0,0,0,0,0,-13.4557338450741,-13.4557338450741,0,0,0,0,0,0,0,0,0,0,-1000114.7537703,0,0,0,0,0,0,0,0,999861.173885226,0,0,0,0,0,0,0,0,0,0,0,999861.173885226,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710335246965512,0,0,0,0,0,0,0,0,0,0,0,-0.137844693702532,0.352285688722168,0,0,-0.0700731140154348,0,0,0,0,0,0,0,-0.0700003629276686,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,537.188544116988,537.188544116988,0,0,0,0,0,0,0,0,0,0,-1004940.22343719,1002006.54362272,0,0,421.587830120486,421.587830120486,0,0,0,0,0,0,633.457524213495,633.457524213495,0,0,0,0,0,0,0,0,0,0,-994601.851407948,0,0,0,0,0,0,0,0,1006530.53920619,0,0,0,0,0,0,0,0,0,0,0,1006530.53920619,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710258099231745,0,0,0,0,0,0,0,0,0,0,0,-0.140178224459363,0.350145679748801,0,0,-0.0700661564728274,0,0,0,0,0,0,0,-0.0700021154913298,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,34.222544931372,34.222544931372,0,0,0,0,0,0,0,0,0,0,-1000314.8020055,993588.567009108,0,0,26.8675595838737,26.8675595838737,0,0,0,0,0,0,40.3435561774616,40.3435561774616,0,0,0,0,0,0,0,0,0,0,-999655.94574482,0,0,0,0,0,0,0,0,1000416.22785372,0,0,0,0,0,0,0,0,0,0,0,1000416.22785372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0709986533422909,0,0,0,0,0,0,0,0,0,0,0,-0.138154590550199,0.347720499406311,0,0,-0.0699607692121066,0,0,0,0,0,0,0,-0.070030813710435,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,6.68975756866653,6.68975756866653,0,0,0,0,0,0,8.95960860656955,8.95960860656955,0,0,-1000046.96374013,997179.219999,0,0,-6.46901364403112,-6.46901364403112,0,0,0,0,0,0,4.10278425519037,4.10278425519037,0,0,0,0,0,0,0,0,0,-999911.615371665,0,0,0,0,0,0,0,0,1000060.24646743,0,0,0,0,0,0,0,0,-999861.173885226,0,0,0,0,0,0,0,0,1000060.24646743,0,0,0,0,0,0,0,0,0,-0.0710362817472667,0,0,0,0,0,0,0,-0.0710730831129618,0,0,0,-0.0680936445282936,0.349008078045094,0,0,-0.0698234010431359,0,0,0,0,0,0,0,-0.0699944953584088,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-314.562059357811,-314.562059357811,0,0,0,0,0,0,-421.19538722176,-421.19538722176,0,0,-997790.968403646,999194.338171343,0,0,304.597421520438,304.597421520438,0,0,0,0,0,0,-192.970297104984,-192.970297104984,0,0,0,0,0,0,0,0,0,-1004157.71399714,0,0,0,0,0,0,0,0,997165.937271698,0,0,0,0,0,0,0,0,-1006530.53920619,0,0,0,0,0,0,0,0,997165.937271698,0,0,0,0,0,0,0,0,0,-0.071038025248108,0,0,0,0,0,0,0,-0.0710676101605021,0,0,0,-0.0670672379316558,0.349936779256152,0,0,-0.0698466303323418,0,0,0,0,0,0,0,-0.0699970640539391,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-20.0577737642986,-20.0577737642986,0,0,0,0,0,0,-26.8659580257674,-26.8659580257674,0,0,-999859.211352149,1002786.54791801,0,0,19.3853267273472,19.3853267273472,0,0,0,0,0,0,-12.2999789623031,-12.2999789623031,0,0,0,0,0,0,0,0,0,-1000264.99440839,0,0,0,0,0,0,0,0,999819.369303291,0,0,0,0,0,0,0,0,-1000416.22785372,0,0,0,0,0,0,0,0,999819.369303291,0,0,0,0,0,0,0,0,0,-0.0709964065867806,0,0,0,0,0,0,0,-0.070960861863767,0,0,0,-0.0677424904912376,0.350989236488403,0,0,-0.0702011111721388,0,0,0,0,0,0,0,-0.0700364502371939,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,18.6430927038996,18.6430927038996,0,0,0,0,0,0,6.46910687871318,6.46910687871318,0,0,-1000160.92860224,990449.959662738,0,0,0,0,0,0,0,0,0,0,17.8895048266275,17.8895048266275,0,0,0,0,0,0,0,0,-999789.737392714,0,0,0,0,0,0,0,0,1000203.92872532,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000060.24646743,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710202818727207,0,0,0,0,0,0,0,-0.0708234019831014,0,0,0,-0.139218987838377,0.346642435838878,0,0,0,0,0,0,0,0,0,0,-0.0700082901082964,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,-876.111886261751,-876.111886261751,0,0,0,0,0,0,-304.390722688526,-304.390722688526,0,0,-992431.703164106,997363.838157946,0,0,0,0,0,0,0,0,0,0,-840.762852169367,-840.762852169367,0,0,0,0,0,0,0,0,-1009890.99351164,0,0,0,0,0,0,0,0,990406.959539653,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-997165.937271698,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710295967498492,0,0,0,0,0,0,0,-0.0708465861152894,0,0,0,-0.135776698301738,0.34978600362384,0,0,0,0,0,0,0,0,0,0,-0.0700089944889218,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,-55.8917431473313,-55.8917431473313,0,0,0,0,0,0,-19.3844920746244,-19.3844920746244,0,0,-999517.502987933,1009438.14766634,0,0,0,0,0,0,0,0,0,0,-53.6308177074652,-53.6308177074652,0,0,0,0,0,0,0,0,-1000630.40843497,0,0,0,0,0,0,0,0,999388.581782399,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999819.369303291,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710120788192635,0,0,0,0,0,0,0,-0.0712011083539383,0,0,0,-0.13855086278332,0.353348474105155,0,0,0,0,0,0,0,0,0,0,-0.0700231926757699,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,13.4561374354576,13.4561374354576,0,0,0,0,0,0,0,0,0,0,-1000138.06501875,992489.833801141,0,0,8.59926580788568,8.59926580788568,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999861.173885226,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1000160.11985362,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710002241015538,0,0,0,0,0,0,0,0,0,0,0,-0.209134137326719,0.347363722638197,0,0,-0.0699217352075871,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,-632.569712929349,-632.569712929349,0,0,0,0,0,0,0,0,0,0,-993506.052050475,998481.656266873,0,0,-404.451809485839,-404.451809485839,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1006530.53920619,0,0,0,0,0,0,0,0,0,0,0,0,0,0,992467.778966278,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.071008646030536,0,0,0,0,0,0,0,0,0,0,0,-0.206419450636998,0.349831975272874,0,0,-0.0699330534272303,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-40.3399440931136,-40.3399440931136,0,0,0,0,0,0,0,0,0,0,-999586.048786897,1007445.70934105,0,0,-25.774404633235,-25.774404633235,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000416.22785372,0,0,0,0,0,0,0,0,0,0,0,0,0,0,999519.92935107,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710312299382888,0,0,0,0,0,0,0,0,0,0,0,-0.208694094893451,0.352629140071907,0,0,-0.0701064864785274,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.10274673552472,-4.10274673552472,0,0,0,0,0,0,-8.59910101799077,-8.59910101799077,0,0,-999974.482308221,1001448.3028934,0,0,7.31738817846815,7.31738817846815,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000060.24646743,0,0,0,0,0,0,0,0,0,0,0,-1000160.11985362,0,0,0,0,0,0,0,0,999969.097526927,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0709945556048762,0,0,0,0,0,0,0,-0.0709215416943785,0,0,0,-0.138893440022622,0.350508790686142,0,0,-0.0701793343212826,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,193.052852618787,193.052852618787,0,0,0,0,0,0,404.815559983562,404.815559983562,0,0,-1001200.43787511,1000345.90145999,0,0,-343.911657271068,-343.911657271068,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-997165.937271698,0,0,0,0,0,0,0,0,0,0,0,-992467.778966278,0,0,0,0,0,0,0,0,1001453.68767469,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0709942299912108,0,0,0,0,0,0,0,-0.0709421565098694,0,0,0,-0.139427508149709,0.350032428081143,0,0,-0.0701595959818201,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12.3003147612105,12.3003147612105,0,0,0,0,0,0,25.7758796498834,25.7758796498834,0,0,-1000076.52678177,998566.372825627,0,0,-21.9484365427487,-21.9484365427487,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999819.369303291,0,0,0,0,0,0,0,0,0,0,0,-999519.92935107,0,0,0,0,0,0,0,0,1000092.65166041,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710362696064972,0,0,0,0,0,0,0,-0.0711070666689879,0,0,0,-0.139031542627119,0.349492586781448,0,0,-0.0698584860884893,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-17.8887914772881,-17.8887914772881,0,0,0,0,0,0,-7.31726876923409,-7.31726876923409,0,0,-999831.693684471,1009077.87574601,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000203.92872532,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999969.097526927,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710084940370217,0,0,0,0,0,0,0,-0.0711795135458062,0,0,0,-0.209744645947434,0.353186078923656,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,842.331535671518,842.331535671518,0,0,0,0,0,0,344.17230700849,344.17230700849,0,0,-1007918.40812958,1001764.86496989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-990406.959539653,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1001453.68767469,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0709994014484615,0,0,0,0,0,0,0,-0.0711511650552634,0,0,0,-0.213017488452914,0.350203066661158,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,53.6372020286134,53.6372020286134,0,0,0,0,0,0,21.9495049627517,21.9495049627517,0,0,-1000504.6112062,990997.141723021,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-999388.581782399,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1000092.65166041,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.0710225812575523,0,0,0,0,0,0,0,-0.070857948738302,0,0,0,-0.210321387896899,0.34682254686407,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1.61232952550381,17.3259782741302,0,0,-0.476042304755734,-0.476042304755734,0,0,0,0,0,0,-0.476529357197454,-0.476529357197454,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.73852394191485,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.73852394191485,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0.0602896370954553,18.9929122524357,0,0,-0.47739318449351,-0.47739318449351,0,0,0,0,0,0,-0.478383920637558,-0.478383920637558,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.72951998552892,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.72951998552892,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-1.57076977684911,17.4954638003806,0,0,-0.474448343813498,-0.474448343813498,0,0,0,0,0,0,-0.474003048546973,-0.474003048546973,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.7623816571778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.7623816571778,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-0.480886374652306,-0.480886374652306,0,0,-1.10079494743586,17.953356602787,0,0,-0.477236262769714,-0.477236262769714,0,0,0,0,0,0,-0.474624490323972,-0.474624490323972,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.76044749564729,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.73852394191485,0,0,0,0,0,0,0,0,4.76044749564729,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-0.479322943072824,-0.479322943072824,0,0,0.487142456369274,19.5456735181573,0,0,-0.478971617697908,-0.478971617697908,0,0,0,0,0,0,-0.473660693010736,-0.473660693010736,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.76745571249605,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.72951998552892,0,0,0,0,0,0,0,0,4.76745571249605,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-0.482971648355281,-0.482971648355281,0,0,-1.13451066656792,17.8211110830998,0,0,-0.47331916866445,-0.47331916866445,0,0,0,0,0,0,-0.475707558877716,-0.475707558877716,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.74207105999484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.7623816571778,0,0,0,0,0,0,0,0,4.74207105999484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-0.482855967074393,-0.482855967074393,0,0,-1.56235325632662,17.5617447854146,0,0,0,0,0,0,0,0,0,0,-0.473247526563853,-0.473247526563853,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.77359980610069,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.76044749564729,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-0.482923055939952,-0.482923055939952,0,0,-0.0267538174078308,19.1077638910622,0,0,0,0,0,0,0,0,0,0,-0.470360090610084,-0.470360090610084,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.79037797916678,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.76745571249605,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-0.481155251543835,-0.481155251543835,0,0,-1.63862129106832,17.2510721105154,0,0,0,0,0,0,0,0,0,0,-0.476944422689672,-0.476944422689672,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.7298818157602,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74207105999484,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-0.48089052889726,-0.48089052889726,0,0,0,0,0,0,0,0,0,0,-1.09152585406242,17.9915336495676,0,0,-0.474628931282051,-0.474628931282051,0,0,0,0,0,0,-0.473664308242654,-0.473664308242654,0,0,0,0,0,0,0,0,0,0,-4.73852394191485,0,0,0,0,0,0,0,0,4.76589094191551,0,0,0,0,0,0,0,0,0,0,0,4.76589094191551,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-0.479129673648487,-0.479129673648487,0,0,0,0,0,0,0,0,0,0,0.476890562701471,19.5668292647523,0,0,-0.473453294484112,-0.473453294484112,0,0,0,0,0,0,-0.471507540526748,-0.471507540526748,0,0,0,0,0,0,0,0,0,0,-4.72951998552892,0,0,0,0,0,0,0,0,4.77692765226587,0,0,0,0,0,0,0,0,0,0,0,4.77692765226587,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-0.482959654428648,-0.482959654428648,0,0,0,0,0,0,0,0,0,0,-1.14435563111168,17.7839411670196,0,0,-0.475694612002936,-0.475694612002936,0,0,0,0,0,0,-0.476576027499387,-0.476576027499387,0,0,0,0,0,0,0,0,0,0,-4.7623816571778,0,0,0,0,0,0,0,0,4.73702662146594,0,0,0,0,0,0,0,0,0,0,0,4.73702662146594,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,-0.48283316844006,-0.48283316844006,0,0,0,0,0,0,-0.483320485766546,-0.483320485766546,0,0,-0.653584426454616,18.3157175915326,0,0,-0.473223509023094,-0.473223509023094,0,0,0,0,0,0,-0.475491862338708,-0.475491862338708,0,0,0,0,0,0,0,0,0,-4.76044749564729,0,0,0,0,0,0,0,0,4.74440686114667,0,0,0,0,0,0,0,0,-4.76589094191551,0,0,0,0,0,0,0,0,4.74440686114667,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,-0.483997921804344,-0.483997921804344,0,0,0,0,0,0,-0.484975560350532,-0.484975560350532,0,0,0.995434940168271,19.9617520683004,0,0,-0.471473203467884,-0.471473203467884,0,0,0,0,0,0,-0.476047300618851,-0.476047300618851,0,0,0,0,0,0,0,0,0,-4.76745571249605,0,0,0,0,0,0,0,0,4.73966640618382,0,0,0,0,0,0,0,0,-4.77692765226587,0,0,0,0,0,0,0,0,4.73966640618382,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,-0.481223454587001,-0.481223454587001,0,0,0,0,0,0,-0.480777271974229,-0.480777271974229,0,0,-0.629513008581732,18.4069504304793,0,0,-0.47701359338184,-0.47701359338184,0,0,0,0,0,0,-0.474941505265021,-0.474941505265021,0,0,0,0,0,0,0,0,0,-4.74207105999484,0,0,0,0,0,0,0,0,4.75693249028916,0,0,0,0,0,0,0,0,-4.73702662146594,0,0,0,0,0,0,0,0,4.75693249028916,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,-0.483996759072047,-0.483996759072047,0,0,0,0,0,0,-0.481386370755409,-0.481386370755409,0,0,-1.16146846434246,17.7262724105024,0,0,0,0,0,0,0,0,0,0,-0.477049897607695,-0.477049897607695,0,0,0,0,0,0,0,0,-4.77359980610069,0,0,0,0,0,0,0,0,4.72889896451955,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74440686114667,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,-0.487020614627862,-0.487020614627862,0,0,0,0,0,0,-0.481625209668949,-0.481625209668949,0,0,0.546596487851598,19.4242349597071,0,0,0,0,0,0,0,0,0,0,-0.479815617656994,-0.479815617656994,0,0,0,0,0,0,0,0,-4.79037797916678,0,0,0,0,0,0,0,0,4.71297138162444,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.73966640618382,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,-0.480188427031373,-0.480188427031373,0,0,0,0,0,0,-0.482569981562914,-0.482569981562914,0,0,-1.08381437949533,18.0312678554781,0,0,0,0,0,0,0,0,0,0,-0.473543714922803,-0.473543714922803,0,0,0,0,0,0,0,0,-4.7298818157602,0,0,0,0,0,0,0,0,4.77129549912168,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75693249028916,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,-0.483311324069199,-0.483311324069199,0,0,0,0,0,0,0,0,0,0,-1.62705601668172,17.285505149804,0,0,-0.475483350274719,-0.475483350274719,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.76589094191551,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.73362747949309,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,-0.485410193009288,-0.485410193009288,0,0,0,0,0,0,0,0,0,0,0.0652596731148528,18.9698883987308,0,0,-0.47644871216546,-0.47644871216546,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.77692765226587,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.72109329610554,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.480805664690484,-0.480805664690484,0,0,0,0,0,0,0,0,0,0,-1.568999026459,17.5220629388667,0,0,-0.474967584057198,-0.474967584057198,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.73702662146594,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.7669165767682,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.481407004842955,-0.481407004842955,0,0,0,0,0,0,-0.480440782561222,-0.480440782561222,0,0,-1.11818380470865,17.9030943926657,0,0,-0.477069671157993,-0.477069671157993,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74440686114667,0,0,0,0,0,0,0,0,0,0,0,-4.73362747949309,0,0,0,0,0,0,0,0,4.75424404839827,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.480657853721154,-0.480657853721154,0,0,0,0,0,0,-0.478784631212474,-0.478784631212474,0,0,0.493801624487004,19.5166109868,0,0,-0.478874387788298,-0.478874387788298,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.73966640618382,0,0,0,0,0,0,0,0,0,0,0,-4.72109329610554,0,0,0,0,0,0,0,0,4.75668708308027,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.482507774882317,-0.482507774882317,0,0,0,0,0,0,-0.483393618315547,-0.483393618315547,0,0,-1.12996586838915,17.8568629555401,0,0,-0.473482065401585,-0.473482065401585,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75693249028916,0,0,0,0,0,0,0,0,0,0,0,-4.7669165767682,0,0,0,0,0,0,0,0,4.74781907458474,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.480034837519482,-0.480034837519482,0,0,0,0,0,0,-0.482302078734292,-0.482302078734292,0,0,-1.57206928715852,17.5423962374812,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.72889896451955,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75424404839827,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.47713222303591,-0.47713222303591,0,0,0,0,0,0,-0.481743846304109,-0.481743846304109,0,0,-0.0291594536385996,19.0948955830024,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.71297138162444,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.75668708308027,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05,0},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.483730274843943,-0.483730274843943,0,0,0,0,0,0,-0.481659493601987,-0.481659493601987,0,0,-1.64224537804411,17.2564892500634,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.77129549912168,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-4.74781907458474,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-8.78778157533803e-05,8.78778157533803e-05},
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}
  };

   AMP::Discretization::DOFManager::shared_ptr left_DOFManager  = testJacobian->getLeftDOFManager();
   AMP::Discretization::DOFManager::shared_ptr right_DOFManager = testJacobian->getRightDOFManager();
   bool equal_to_leftDOFManager, equal_to_rightDOFManager;
   if (*left_DOFManager  == *subchannelDOFManager) equal_to_leftDOFManager  = true;
   if (*right_DOFManager == *subchannelDOFManager) equal_to_rightDOFManager = true;
   AMP_ASSERT(equal_to_leftDOFManager);
   AMP_ASSERT(equal_to_rightDOFManager);

   // create global ID maps
   std::map<size_t,AMP::Mesh::MeshElement> elements_by_globalID;
   std::map<size_t,size_t>                 variables_by_globalID;
   createGlobalIDMaps(subchannelDOFManager,subchannelMesh,elements_by_globalID,variables_by_globalID);
   
double max_entry = 0.0;
for (size_t i = 0; i < num_dofs_AMP; ++i)
   for (size_t j = 0; j < num_dofs_AMP; ++j) {
      double entry = testJacobian->getValueByGlobalID(i,j);
      max_entry = std::max(max_entry,entry);
   }
std::cout << "max value: " << max_entry << std::endl;
  // initialize success boolean for known residual comparison test
  /*
  bool passed_known_test = JacobianIsCorrect(testJacobian,
                                              knownJacobian,
                                              subchannelDOFManager,
                                              subchannelMesh,
                                              interiorLateralFaceMap,
                                              elements_by_globalID,
                                              variables_by_globalID);

  if (passed_known_test) ut->passes(exeName+": known value test");
  else ut->failure(exeName+": known residual test");
*/
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
        "testSubchannelFourEqLinearOperator"
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


