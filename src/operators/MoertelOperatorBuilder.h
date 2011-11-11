#ifndef included_MoertelOperatorBuilder_h
#define included_MoertelOperatorBuilder_h

#include <map>

#include "EpetraMatrixOperator.h"
#include "MoertelOperatorBuilderParameters.h"

#include "mrtr_segment_bilinearquad.H"
#include "mrtr_manager.H"

#include "ampmesh/Mesh.h"
#include "descritization/DOF_Manager.h"


namespace AMP {
namespace Operator {
/*
  class MoertelOperatorBuilder
  {
    private:
      Epetra_MpiComm                  d_Comm;
      MOERTEL::Manager                d_Manager;
      Epetra_CrsMatrix               *d_M;
      Epetra_CrsMatrix               *d_D;
      Epetra_CrsMatrix               *d_MT;
      Epetra_CrsMatrix               *d_DT;
      AMP::LinearAlgebra::Vector::shared_ptr              d_LambdaVector;
      size_t                          d_NumNodes;
      boost::shared_ptr<Database>     d_DB;
      std::map<int,int>               d_NodeRenumber;
      int                             d_CurLocalNodeID;

      typedef AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator   NodeIterator;
      typedef AMP::Mesh::MeshManager::Adapter::BoundarySideIterator        SideIterator;

      void  subMatrix ( const Epetra_CrsMatrix *in , Epetra_CrsMatrix *&out , int row_start , int row_end , int col_start , int col_end );
      void  transposeMatrix ( const Epetra_CrsMatrix * , Epetra_CrsMatrix *& );
      void  getNodes ( AMP::Mesh::Mesh::shared_ptr , NodeIterator , NodeIterator , MOERTEL::Interface & , AMP::Discretization::DOFManager::shared_ptr , int , size_t );

      size_t buildInformation ( AMP::Mesh::MeshManager::Adapter::shared_ptr , short int , AMP::Mesh::DOFMap::shared_ptr , MOERTEL::Interface & , int , size_t , size_t );
      void  buildLambdaVector ();


    public:
      MoertelOperatorBuilder ( boost::shared_ptr<MoertelOperatorBuilderParameters> );
      ~MoertelOperatorBuilder ();

      boost::shared_ptr<Operator>   createMOperator ( AMP::LinearAlgebra::Variable::shared_ptr , AMP::LinearAlgebra::Variable::shared_ptr );
      boost::shared_ptr<Operator>   createDOperator ( AMP::LinearAlgebra::Variable::shared_ptr , AMP::LinearAlgebra::Variable::shared_ptr );
      boost::shared_ptr<Operator>   createMTOperator ( AMP::LinearAlgebra::Variable::shared_ptr , AMP::LinearAlgebra::Variable::shared_ptr );
      boost::shared_ptr<Operator>   createDTOperator ( AMP::LinearAlgebra::Variable::shared_ptr , AMP::LinearAlgebra::Variable::shared_ptr );

      AMP::LinearAlgebra::Vector::shared_ptr            createLambdaVector ( const std::string &name = "" )
      {
        return d_LambdaVector->cloneVector ( name );
      }
  };

*/

}
}

#endif
