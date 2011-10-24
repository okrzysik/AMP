#ifndef included_AMP_VectorVariable_h
#define included_AMP_VectorVariable_h

#include "Variable.h"

namespace AMP {
namespace LinearAlgebra {

  template <typename VAR_TYPE , int N>
    class VectorVariable : public VAR_TYPE
  {
    public:
      VectorVariable ( std::string name , typename VAR_TYPE::Discretization::shared_ptr p = typename VAR_TYPE::Discretization::shared_ptr () ) : VAR_TYPE ( name , p ) {}
      virtual  size_t   variableID () const { return N + VAR_TYPE::OFFSET; }
      virtual  Variable::shared_ptr   cloneVariable ( const std::string &name ) const
      { 
        VectorVariable<VAR_TYPE,N> *retVal = new VectorVariable<VAR_TYPE,N> ( name ); 
        retVal->setDiscretization ( *this );
        retVal->setUnits ( this->getUnits() );
        return Variable::shared_ptr ( retVal );
      }

      virtual bool  operator == ( const Variable &rhs ) const
      {
        if ( rhs.isA<VectorVariable<VAR_TYPE,N> >() )
          return VAR_TYPE::operator == ( rhs );
        return false;
      }

      virtual  size_t   DOFsPerObject () const { return N; }
  };

}
}

#endif
