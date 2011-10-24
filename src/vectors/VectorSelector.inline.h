namespace AMP {
namespace LinearAlgebra {

  inline
  VectorSelector::~VectorSelector ()
  {
  }

  inline
  bool  VectorSelector::isSelected ( Vector::const_shared_ptr ) const
  {
    return true;
  }

  inline
  VectorSubsetter::shared_ptr  VectorSelector::getSubsetter () const
  {
    return VectorSubsetter::shared_ptr ( new VectorSubsetter () );
  }
  
  inline
  VS_ByVariableName::~VS_ByVariableName ()
  {
  }

  inline
  VS_ByVariableName::VS_ByVariableName ( std::string  n )
       : d_VecName ( n ) 
  {
  }

  inline
  bool   VS_ByVariableName::isSelected ( Vector::const_shared_ptr v ) const
  {
    return v->getVariable()->getName() == d_VecName;
  }

  inline
  VS_Stride::VS_Stride ( const std::string &n , size_t a , size_t b ) : 
    d_Offset ( a ) , 
    d_Stride ( b ) , 
    d_Name ( n ) 
  {
  }

  inline
  VS_Stride::~VS_Stride () 
  {
  }

  inline
  VectorSubsetter::shared_ptr  VS_Stride::getSubsetter () const
  {
    return VectorSubsetter::shared_ptr ( new VectorStriderSubsetter ( d_Name , d_Offset , d_Stride ) );
  }

}
}

