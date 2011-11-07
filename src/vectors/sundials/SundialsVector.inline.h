
namespace AMP {
namespace LinearAlgebra {


  inline
  SundialsVector::SundialsVector ( ) 
  {
  }


  inline
	N_Vector &SundialsVector::getNVector()            
  { 
    return d_n_vector; 
  }

  inline
	const N_Vector &SundialsVector::getNVector() const
  { 
    return d_n_vector; 
  }

}
}


