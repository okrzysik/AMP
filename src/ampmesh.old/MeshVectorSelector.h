#ifndef  included_AMP_MeshVectorSelector
#define  included_AMP_MeshVectorSelector

#include "vectors/VectorSelector.h"

#include "MeshVariable.h"

namespace AMP { 
namespace Mesh {

  template <typename T>
  class VS_ByMeshNameTmpl : public AMP::LinearAlgebra::VectorSelector
  {
    protected:
      std::string   d_MeshName;

    public:
      virtual ~VS_ByMeshNameTmpl () {}

      VS_ByMeshNameTmpl ( std::string n ) : d_MeshName ( n ) {}

      virtual bool  isSelected ( AMP::LinearAlgebra::Vector::const_shared_ptr v ) const
      {
        if ( v->getVariable()->isA<MeshVariable <T> >() )
          if ( v->getVariable()->castTo<MeshVariable <T> >().getMesh() )
            if ( v->getVariable()->castTo<MeshVariable <T> >().getMesh()->getMeshName() == d_MeshName )
              return true;
        return false;
      }
  };

  template <typename T>
  class VS_ByMeshTmpl : public AMP::LinearAlgebra::VectorSelector
  {
    protected:
      typename T::shared_ptr   d_Mesh;

    public:
      virtual ~VS_ByMeshTmpl () {}

      VS_ByMeshTmpl ( typename T::shared_ptr  p ) : d_Mesh ( p )  {}

      virtual bool  isSelected ( AMP::LinearAlgebra::Vector::const_shared_ptr v ) const
      {
        if ( v->getVariable()->isA<MeshVariable <T> >() )
          if ( v->getVariable()->castTo<MeshVariable <T> >().getMesh() )
            if ( v->getVariable()->castTo<MeshVariable <T> >().getMesh() == d_Mesh.get() )
              return true;
        return false;
      }
  };

  template <typename MESH,typename ITER>
  class VS_ByMeshIteratorTmpl : public VS_ByMeshTmpl<MESH>
  {
    protected:
      typename  MESH::shared_ptr  d_Mesh;
                ITER              d_Begin;
                ITER              d_End;
                size_t            d_DOFsPerObj;

    public:
      VS_ByMeshIteratorTmpl ( typename MESH::shared_ptr  m , ITER b , ITER e , size_t dpo )
          :  VS_ByMeshTmpl<MESH> ( m ) , 
             d_Begin ( b ) , 
             d_End ( e ) , 
             d_DOFsPerObj ( dpo )
      {
      }

      virtual  AMP::LinearAlgebra::VectorSubsetter::shared_ptr  getSubsetter () const;
  };


  template <typename MESH , typename ITER>
  AMP::LinearAlgebra::VectorSubsetter::shared_ptr   VS_ByMeshIteratorTmpl<MESH,ITER>::getSubsetter () const
  {
    AMP::LinearAlgebra::VectorRandomAccessSubsetter *retval = new AMP::LinearAlgebra::VectorRandomAccessSubsetter ( "temp" );
    ITER t = d_Begin;
    while ( t != d_End )
    {
      for ( size_t i = 0 ; i != d_DOFsPerObj ; i++ )
      {
        retval->addID ( t->globalID() * d_DOFsPerObj + i );
      }
      t++;
    }
    return AMP::LinearAlgebra::VectorSubsetter::shared_ptr ( retval );
  }

}
}


#endif
