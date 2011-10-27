namespace AMP { 
namespace Mesh {

  inline
  LibMeshPoint::LibMeshPoint ( ::Point &rhs ) : d_Pt ( rhs ) 
  {
  }

  inline
  LibMeshPoint::~LibMeshPoint () 
  {
  }

  inline
  double    LibMeshPoint::x () const 
  { 
    return d_Pt(0); 
  }

  inline
  double    LibMeshPoint::y () const 
  { 
    return d_Pt(1); 
  }

  inline
  double    LibMeshPoint::z () const 
  { 
    return d_Pt(2); 
  }

  inline
  double   &LibMeshPoint::x ()  
  { 
    return d_Pt(0); 
  }

  inline
  double   &LibMeshPoint::y ()  
  { 
    return d_Pt(1); 
  }

  inline
  double   &LibMeshPoint::z ()  
  { 
    return d_Pt(2); 
  }

  inline
  double  &LibMeshPoint::operator ()( size_t i )       
  { 
    return d_Pt(i); 
  }

  inline
  double   LibMeshPoint::operator ()( size_t i ) const 
  { 
    return d_Pt(i); 
  }

}
}

