#ifndef included_AMP_MeshVariable_h
#define included_AMP_MeshVariable_h

#include "vectors/VectorVariable.h"
#include <boost/shared_ptr.hpp>
#include "MeshType.h"

namespace AMP { 
namespace Mesh {

  // These are stubs for use when we incorporate more operators.

  class MeshDOFObject : public Castable
  {
    public:
      typedef boost::shared_ptr<MeshDOFObject>   shared_ptr;

      MeshDOFObject () {}
      virtual ~MeshDOFObject () {}

      virtual  size_t   numDOFPoints () = 0;
  };

  template <typename MESHTYPE>
  class MeshVariable : public AMP::LinearAlgebra::Variable
  {
    public:
      typedef MESHTYPE   Discretization ;
    protected:
      // This should be a shared_ptr, except it can introduce a cyclic dependency.
      MESHTYPE     *d_Mesh;

      void setDiscretization ( const MeshVariable <MESHTYPE> &rhs )
      {
        d_Mesh = rhs.d_Mesh; 
      }

    public:
      virtual ~MeshVariable() {}
      MeshVariable ( std::string name , typename MESHTYPE::shared_ptr p = typename MESHTYPE::shared_ptr() ) : Variable ( name ) , d_Mesh ( p.get() ) {}

      void  setMesh ( typename MESHTYPE::shared_ptr m )
      {
        d_Mesh = m.get();
      }

      MESHTYPE  *getMesh ()
      {
        return d_Mesh;
      }

      virtual bool  operator == ( const Variable &rhs ) const
      {
        bool retVal = false;
        if ( rhs.isA<MeshVariable<MESHTYPE> >() )
        {
          const MeshVariable<MESHTYPE>  &rhs_t = rhs.castTo<MeshVariable<MESHTYPE> >();
          retVal = d_VariableName == rhs_t.d_VariableName;
          retVal &= (d_Mesh == rhs_t.d_Mesh);
        }
        return retVal;
      }

  };

  class NodalVariable : public MeshVariable<MESH_TYPE>
  {
    public:
      enum { OFFSET = 0 };
      NodalVariable ( std::string name , MESH_TYPE::shared_ptr p = MESH_TYPE::shared_ptr()) : MeshVariable<MESH_TYPE> ( name , p ) {}
      virtual bool operator == ( const Variable &rhs ) const
      {
        return MeshVariable<MESH_TYPE>::operator == ( rhs );
      }
  };

  class IntegrationPointVariable : public MeshVariable<MESH_TYPE>
  {
    public:
      enum { OFFSET = 1024 };
      IntegrationPointVariable ( std::string name , MESH_TYPE::shared_ptr p = MESH_TYPE::shared_ptr()) : MeshVariable<MESH_TYPE> ( name , p ) {}
      virtual bool operator == ( const Variable &rhs ) const
      {
        return MeshVariable<MESH_TYPE>::operator == ( rhs );
      }
  };


  // Popular variable types:
  typedef AMP::LinearAlgebra::VectorVariable<NodalVariable, 1>               NodalScalarVariable;
  typedef AMP::LinearAlgebra::VectorVariable<NodalVariable, 2>               Nodal2VectorVariable;
  typedef AMP::LinearAlgebra::VectorVariable<NodalVariable, 3>               Nodal3VectorVariable;

  typedef AMP::LinearAlgebra::VectorVariable<IntegrationPointVariable, 1>    SingleGaussPointVariable;

  class RunTimeIntegrationPointVariable : public IntegrationPointVariable
  {
    private:
      size_t    d_NumEntriesPerObject;

    public:

      RunTimeIntegrationPointVariable ( 
                 std::string name , 
                 size_t NumEntriesPerObject , 
                 MESH_TYPE::shared_ptr p = MESH_TYPE::shared_ptr() ) 
        : IntegrationPointVariable ( name , p ) , d_NumEntriesPerObject ( NumEntriesPerObject ) {}
      virtual size_t   variableID() const  { return OFFSET + d_NumEntriesPerObject; }
      virtual  AMP::LinearAlgebra::Variable::shared_ptr   cloneVariable ( const std::string &name ) const
      {
        AMP::LinearAlgebra::Variable::shared_ptr retVal ( new RunTimeIntegrationPointVariable ( name , d_NumEntriesPerObject ) );
        return retVal;
      }

      virtual bool   operator == ( const Variable &rhs ) const
      {
        if ( !rhs.isA<IntegrationPointVariable> () )
          return false;
        if ( !IntegrationPointVariable::operator == ( rhs ) )
          return false;
        return d_NumEntriesPerObject == rhs.castTo<RunTimeIntegrationPointVariable>().d_NumEntriesPerObject;
      }

      virtual size_t  DOFsPerObject () const { return d_NumEntriesPerObject; }
  };

}
}

#endif
