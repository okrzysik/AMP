namespace AMP {
namespace LinearAlgebra {

  inline
  VectorSubsetter::~VectorSubsetter () 
  {
  }

  inline
  Vector::shared_ptr  VectorSubsetter::subset ( Vector::shared_ptr p ) 
  { 
    return p; 
  }


  inline
  VectorStriderSubsetter::VectorStriderSubsetter ( const std::string &n , size_t a , size_t b )
   : d_Offset ( a ) , d_Stride ( b ) , d_Name ( n ) 
  {
  }

  inline
  VectorStriderSubsetter::~VectorStriderSubsetter () 
  {
  }

  inline
  Vector::shared_ptr  VectorStriderSubsetter::subset ( Vector::shared_ptr p )
  { 
    Variable::shared_ptr  variable ( new StridedVariable( d_Name , d_Offset , d_Stride ) );
    Vector::shared_ptr  vector = SubsetVector::view ( p, variable ); 
    return vector;
  }

}
}

