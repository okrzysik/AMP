#include "operators/map/Map3Dto1D.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

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

    AMP_INSIST( myparams->d_db->keyExists("InputVariable"), "key not found" );
    std::string inpVar = myparams->d_db->getString("InputVariable");
    d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar) );

    AMP_INSIST( myparams->d_db->keyExists("OutputVariable"), "key not found" );
    std::string outVar = myparams->d_db->getString("OutputVariable");
    d_outVariable.reset(new AMP::LinearAlgebra::Variable(outVar) );

 } 

 void Map3Dto1D :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
     AMP::LinearAlgebra::Vector::shared_ptr  &r, const double , const double )
 { 
   AMP_ASSERT(u != NULL);

   AMP::LinearAlgebra::Vector::shared_ptr inputVec = u->subsetVectorForVariable(d_inpVariable);
   inputVec->makeConsistent ( AMP::LinearAlgebra::Vector::CONSISTENT_SET );
   AMP_ASSERT( inputVec != NULL);

   if(d_iDebugPrintInfoLevel>5)
   {
     AMP::pout << "The input to Map3Dto1D " << std::endl;
     AMP::pout << inputVec << std::endl;
   }

//   AMP::LinearAlgebra::Vector::shared_ptr outputVec =  r->subsetVectorForVariable(d_outVariable);

   AMP_ASSERT(outputVec  != NULL);

   const unsigned int numPoints = outputVec->getLocalSize();

AMP_ERROR("Not finished");
/*
   AMP::Mesh::MeshIterator bnd = d_Mesh->beginSideBoundary( d_boundaryId );
   AMP::Mesh::MeshIterator end_bnd = bnd.end();
   //AMP::Mesh::MeshManager::Adapter::BoundarySideIterator bnd = d_MeshAdapter->beginSideBoundary( d_boundaryId );
   //AMP::Mesh::MeshManager::Adapter::BoundarySideIterator end_bnd = d_MeshAdapter->endSideBoundary( d_boundaryId );

   std::vector<double> mapValues(numPoints,0);
   std::vector<int> numFaceNodes(numPoints,0);

   // Iterator for the solid-clad boundary
   for( ; bnd != end_bnd; ++bnd) {

     AMP::Mesh::MeshElement cur_side = *bnd;
     const unsigned int num_nodes = cur_side.numNodes();

     std::vector<double> zcoords;

     for(unsigned int i=0; i <num_nodes; i++ ) {
       zcoords.push_back(d_MeshAdapter->getNode(cur_side.getNodeID(i)).z() );
     }

     std::sort(zcoords.begin(),zcoords.end());

     std::vector<double> tmpZcoords = zcoords;
     std::vector<int>   tmpIds(zcoords.size());

     for(unsigned int i = 0; i < num_nodes; i++) {
       tmpIds[i] = i;
     }

     std::vector<int> originalNodeOrder(zcoords.size());

     for(unsigned int i=0; i <num_nodes; i++ ) {
       double myZ = d_MeshAdapter->getNode(cur_side.getNodeID(i)).z() ;
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

     z[0] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[0])).z() ;
     z[2] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[2])).z() ;
     z[3] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[3])).z() ;

     y[0] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[0])).y() ;
     y[2] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[2])).y() ;
     y[3] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[3])).y() ;

     x[0] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[0])).x() ;
     x[2] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[2])).x() ;
     x[3] =  d_MeshAdapter->getNode(cur_side.getNodeID(originalNodeOrder[3])).x() ;

     int pickId;     
     if( pow((y[0]-y[3]),2)+pow((x[0]-x[3]),2) < pow((y[0]-y[2]),2)+pow((x[0]-x[2]),2))
     {
       pickId = 3;
     } else {
       pickId = 2;
     }

     // Iterator for the fluid boundary
     for(unsigned int i = 0 ; i < numPoints; i++) {

       double cur_node  = d_zLocations[i];

       // Section of the Clad boundary map corresponding to the fluid Element 
       if( cur_node >= z[0] && cur_node <= z[pickId])
       {
         std::vector<unsigned int> bndGlobalIds;
         dof_map->getDOFs(cur_side, bndGlobalIds);

         mapValues[i]  += ((inputVec)->getValueByGlobalID(bndGlobalIds[originalNodeOrder[0]])* (z[pickId]-cur_node)+ (inputVec)->getValueByGlobalID(bndGlobalIds[originalNodeOrder[pickId]])* (cur_node-z[0]))/(z[pickId]-z[0]);  

         numFaceNodes[i] += 1;

       }

     }//end for i
   }//end for bnd

   std::vector<double>  aggMapValues ( numPoints );
   std::vector<int>  aggNumFaceNodes ( numPoints );
   AMP_MPI myComm(d_MeshAdapter->getComm());
   myComm.sumReduce( (double*) &(mapValues[0]), (double*) &(aggMapValues[0]), numPoints );
   myComm.sumReduce( (int*) &(numFaceNodes[0]), (int*) &(aggNumFaceNodes[0]), numPoints );

   for(unsigned int i = 0 ; i < numPoints; i++) {
     outputVec->setValueByLocalID ( i , aggMapValues[i]/static_cast<double>(aggNumFaceNodes[i]) );
   }

   if(d_iDebugPrintInfoLevel>4)
   {
     AMP::pout << "The output to Map3Dto1D " << std::endl;
     AMP::pout << outputVec << std::endl;
   }
*/

}


}
}

