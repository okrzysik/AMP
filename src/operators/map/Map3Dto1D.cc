
#include "operators/map/Map3Dto1D.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"
#include "discretization/DOF_Manager.h"

namespace AMP {
namespace Operator {


    Map3Dto1D :: Map3Dto1D(const boost::shared_ptr<OperatorParameters>& params): 
      MapOperator (params)
    {
      boost::shared_ptr<MapOperatorParameters> myparams = 
        boost::dynamic_pointer_cast<MapOperatorParameters>(params);
      reset(myparams);
    }


    void Map3Dto1D :: reset(const boost::shared_ptr<OperatorParameters>& params)
    {
      boost::shared_ptr<MapOperatorParameters> myparams =
        boost::dynamic_pointer_cast<MapOperatorParameters>(params);

      AMP_INSIST( ((myparams.get()) != NULL), "NULL parameter" );
      AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );
      AMP_INSIST( !myparams->d_MapComm.isNull(), "NULL communicator" );
      d_Mesh = myparams->d_Mesh;
      d_MapMesh = myparams->d_MapMesh;
      d_MapComm = myparams->d_MapComm;
      AMP_INSIST( d_MapComm.sumReduce<int>(d_MapMesh.get()!=NULL?1:0)>0, "Somebody must own the mesh");

      AMP_INSIST( myparams->d_db->keyExists("InputVariable"), "key not found" );
      std::string inpVar = myparams->d_db->getString("InputVariable");
      d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar) );

      AMP_INSIST( myparams->d_db->keyExists("OutputVariable"), "key not found" );
      std::string outVar = myparams->d_db->getString("OutputVariable");
      d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar) );

    } 


    void Map3Dto1D :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
        AMP::LinearAlgebra::Vector::shared_ptr &r, const double , const double )
    { 

      const unsigned int numPoints = outputVec->getLocalSize();
      std::vector<double> mapValues(numPoints,0);
      std::vector<int> numFaceNodes(numPoints,0);

      // Get the local contributions to the map
      if ( d_MapMesh != NULL ) {
        AMP_ASSERT(u != NULL);
        AMP::LinearAlgebra::Vector::shared_ptr inputVec = subsetInputVector( u );
        inputVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
        AMP_ASSERT( inputVec != NULL);

        AMP::Discretization::DOFManager::shared_ptr dof_map = inputVec->getDOFManager();

        if(d_iDebugPrintInfoLevel>5) {
          AMP::pout << "The input to Map3Dto1D " << std::endl;
          AMP::pout << inputVec << std::endl;
        }

        // AMP::LinearAlgebra::Vector::shared_ptr outputVec =  u->subsetVectorForVariable(d_outpVariable);    // Output vector is a simple vector

        AMP_ASSERT(outputVec  != NULL);

        // Get an iterator over the side elements
        AMP::Mesh::MeshIterator bnd = d_MapMesh->getBoundaryIDIterator( AMP::Mesh::Face, d_boundaryId, 0 );
        AMP::Mesh::MeshIterator end_bnd = bnd.end();

        // Iterator for the solid-clad boundary
        for( ; bnd != end_bnd; ++bnd) {

          AMP::Mesh::MeshElement cur_side = *bnd;
          std::vector<AMP::Mesh::MeshElement> nodes = cur_side.getElements(AMP::Mesh::Vertex);
          AMP_ASSERT(nodes.size()==4);

          std::vector<double> zcoords;
          for(size_t i=0; i <nodes.size(); i++) {
            std::vector<double> coord = nodes[i].coord();
            zcoords.push_back(coord[2]);
          }

          std::sort(zcoords.begin(),zcoords.end());

          std::vector<double> tmpZcoords = zcoords;
          std::vector<int>   tmpIds(zcoords.size());

          for(size_t i=0; i<nodes.size(); i++) {
            tmpIds[i] = i;
          }

          std::vector<int> originalNodeOrder(zcoords.size());

          for(size_t i=0; i<nodes.size(); i++) {
            std::vector<double> coord = nodes[i].coord();
            double myZ = coord[2];
            for(unsigned int j=0; j < tmpZcoords.size(); j++ ) {
              if( fabs(tmpZcoords[j]-myZ) <= 1.e-12 ) {
                originalNodeOrder[tmpIds[j]] = i;
                tmpZcoords.erase(tmpZcoords.begin() + j);  
                tmpIds.erase(tmpIds.begin() + j);
                break;
              }
            }
          }

          std::vector<double> z(4,0);
          std::vector<double> y(4,0);
          std::vector<double> x(4,0);
          for (int i=0; i<4; i++) {
            std::vector<double> coord = nodes[originalNodeOrder[i]].coord();
            x[i] = coord[0];
            y[i] = coord[1];
            z[i] = coord[2];
          }

          int pickId;
          if( pow((y[0]-y[3]),2)+pow((x[0]-x[3]),2) < pow((y[0]-y[2]),2)+pow((x[0]-x[2]),2)) {
            pickId = 3;
          } else {
            pickId = 2;
          }

          // Iterator for the fluid boundary
          for(unsigned int i = 0 ; i < numPoints; i++) {

            double cur_node  = d_zLocations[i];

            // Section of the Clad boundary map corresponding to the fluid Element 
            if( cur_node >= z[0] && cur_node <= z[pickId]) {
              std::vector<size_t> dof1;
              std::vector<size_t> dof2;
              dof_map->getDOFs( nodes[originalNodeOrder[0]].globalID(), dof1 );
              dof_map->getDOFs( nodes[originalNodeOrder[pickId]].globalID(), dof2 );
              AMP_ASSERT(dof1.size()==1&&dof2.size()==1);

              mapValues[i]  += ((inputVec)->getValueByGlobalID(dof1[0]) * (z[pickId]-cur_node) + (inputVec)->getValueByGlobalID(dof2[0])* (cur_node-z[0]))/(z[pickId]-z[0]);  

              numFaceNodes[i] += 1;

            }

          }//end for i
        }//end for bnd
      }

      // Gather the results from all processors
      std::vector<double>  aggMapValues ( numPoints );
      std::vector<int>  aggNumFaceNodes ( numPoints );
      d_MapComm.sumReduce( (double*) &(mapValues[0]), (double*) &(aggMapValues[0]), numPoints );
      d_MapComm.sumReduce( (int*) &(numFaceNodes[0]), (int*) &(aggNumFaceNodes[0]), numPoints );

      // Store the results
      for(unsigned int i = 0 ; i < numPoints; i++) {
        outputVec->setValueByLocalID ( i , aggMapValues[i]/static_cast<double>(aggNumFaceNodes[i]) );
      }

      if(d_iDebugPrintInfoLevel>4) {
        AMP::pout << "The output to Map3Dto1D " << std::endl;
        AMP::pout << outputVec << std::endl;
      }

    }


}
}

