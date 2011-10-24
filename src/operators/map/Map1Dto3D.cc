#include "utils/AMP_MPI.h"
#include "operators/map/Map1Dto3D.h"
#include "ampmesh/MeshUtils.h"
#include "utils/Utilities.h"
#include "utils/InputDatabase.h"

namespace AMP {
namespace Operator {

 void Map1Dto3D :: reset(const boost::shared_ptr<OperatorParameters>& params)
 {
    boost::shared_ptr<MapOperatorParameters> myparams =
         boost::dynamic_pointer_cast<MapOperatorParameters>(params);

    AMP_INSIST( ((myparams.get()) != NULL), "NULL parameter" );
    AMP_INSIST( (((myparams->d_db).get()) != NULL), "NULL database" );

    computeZLocations();

    AMP_INSIST( myparams->d_db->keyExists("InputVariable"), "key not found" );
    std::string inpVar = myparams->d_db->getString("InputVariable");
    d_inpVariable.reset(new AMP::LinearAlgebra::Variable(inpVar));

    AMP_INSIST( myparams->d_db->keyExists("OutputVariable"), "key not found" );
    std::string outVar = myparams->d_db->getString("OutputVariable");
    d_outVariable.reset(new AMP::Mesh::NodalScalarVariable(outVar,d_MapAdapter));

 }

 double *deRefVectorWithoutCompileError ( std::vector<double> &v )
 {
   if ( v.size() == 0 ) return 0;
   return &(v[0]);
 }

 void Map1Dto3D::computeZLocations(){

   AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator bnd = d_MapAdapter->beginBoundary( d_boundaryId );
   AMP::Mesh::MeshManager::Adapter::BoundaryNodeIterator end_bnd = d_MapAdapter->endBoundary( d_boundaryId );

   d_zLocations.clear();
   double Xx=0;
   double Yy=0;
   std::vector<double> t_zLocations;
   if(bnd!=end_bnd)
   {
     t_zLocations.push_back(bnd->z());
     Xx = bnd->x();
     Yy = bnd->y();
     bnd++;
   }

   for( ; bnd != end_bnd; ++bnd) {
     if( (fabs(Xx - (bnd->x())) <= 1.e-12) && (fabs(Yy - (bnd->y())) <= 1.e-12) ){ 
       t_zLocations.push_back(bnd->z());
     }
   }

   //Make Z locations consistent across all processors.
   unsigned int myLen , totLen;
   myLen = t_zLocations.size();
   totLen = d_MapAdapter->getComm().sumReduce(myLen);
   d_zLocations.resize ( totLen );

   // This will crash in debug mode if d_zLocations has no data
   AMP_MPI comm = AMP_MPI( d_MapAdapter->getComm() );
   comm.allGather ( deRefVectorWithoutCompileError ( t_zLocations ) , myLen , &(d_zLocations[0]) );

   // This will create a list of unique z coordinates
   std::sort(d_zLocations.begin(),d_zLocations.end());
   std::vector<double>::iterator  curLook , curEnd;
   curLook = d_zLocations.begin();
   curEnd = d_zLocations.end();
   curEnd--;
   while ( curLook != curEnd )
   {
     std::vector<double>::iterator peek = curLook;
     peek++;
     if ( fabs ( *peek - *curLook ) < 1.e-12 )
     {
       while ( peek != curEnd )
       {
         std::swap ( *peek , *(peek + 1) );
         peek++;
       }
       curEnd--;
     }
     else
     {
       curLook++;
     }
   }
   size_t newSize = curEnd - d_zLocations.begin() + 1;
   d_zLocations.resize ( newSize );

 }

 void Map1Dto3D :: apply(const AMP::LinearAlgebra::Vector::shared_ptr &, const AMP::LinearAlgebra::Vector::shared_ptr &u,
     AMP::LinearAlgebra::Vector::shared_ptr  &r, const double , const double )
 { 
   AMP::Mesh::DOFMap::shared_ptr dof_map = d_MapAdapter->getDOFMap(d_outVariable);

   AMP_ASSERT(u != NULL);

   AMP::LinearAlgebra::Vector::shared_ptr inputVec = u->subsetVectorForVariable(d_inpVariable);

   AMP_ASSERT(inputVec != NULL);

//   AMP::LinearAlgebra::Vector::shared_ptr outputVec =  r->subsetVectorForVariable(d_outVariable);

   AMP_ASSERT(outputVec != NULL);

   //     outputVec->zero();

   std::vector<int> numFaceNodes(outputVec->getLocalSize(),0);

   const unsigned int numPoints = inputVec->getLocalSize();

   std::vector<unsigned int> dofs;
   dofs.push_back(0);

   // Iterator for the fluid boundary
   for(unsigned int i = 0 ; i < numPoints; i++) {


     AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd     = d_MapAdapter->beginOwnedBoundary( d_boundaryId );
     AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd = d_MapAdapter->endOwnedBoundary( d_boundaryId );

     // Iterator for the solid-clad boundary
     for( ; bnd != end_bnd; ++bnd) {

       std::vector<unsigned int> bndGlobalIds;
       dof_map->getDOFs(*bnd, bndGlobalIds, dofs);
       if(fabs(d_zLocations[i]-bnd->z()) <= 1.e-12)
       {
         outputVec->setValueByGlobalID(bndGlobalIds[0], inputVec->getValueByLocalID(i) );
       }

     }//end for bnd
   }//end for i

   if(d_iDebugPrintInfoLevel>4)
   {
     AMP::pout << "The input to Map1Dto3D " << std::endl;
     AMP::pout << inputVec << std::endl;
   }

   outputVec->makeConsistent( AMP::LinearAlgebra::Vector::CONSISTENT_SET );

   if(d_iDebugPrintInfoLevel>5)
   {
     AMP::pout << "The output to Map1Dto3D " << std::endl;
     AMP::pout << outputVec << std::endl;
   }

 }

}
}

