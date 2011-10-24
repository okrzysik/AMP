
#include <stdexcept>
#include <algorithm>
#include <math.h>

#include "utils/Utilities.h"

#include "MultiVector.h"

namespace AMP {
namespace LinearAlgebra {

  inline
  Vector::iterator    MultiVector::begin() 
  { 
    return Vector::iterator ( d_vVectors.front()->begin() , this , 0 ); 
  }

  inline
  Vector::iterator    MultiVector::end() 
  { 
    return Vector::iterator ( d_vVectors.back()->end() , this , d_vVectors.size()-1 ); 
  }

  inline
  Vector::const_iterator    MultiVector::begin() const 
  { 
    return Vector::const_iterator ( d_vVectors.front()->castTo<const Vector>().begin() , this , 0 ); 
  }

  inline
  Vector::const_iterator    MultiVector::end() const 
  { 
    return Vector::const_iterator ( d_vVectors.back()->castTo<const Vector>().end() , this , d_vVectors.size()-1 ); 
  }

  inline
  Vector::iterator    MultiVector::getIterator ( size_t i )
  {
    if ( i < d_vVectors.size() )
      return d_vVectors[i]->begin();
    else
      return d_vVectors.back()->end();
  }

  inline
  Vector::const_iterator    MultiVector::getIterator ( size_t i ) const
  {
    if ( i < d_vVectors.size() )
      return boost::dynamic_pointer_cast<const Vector> (d_vVectors[i])->begin();
    else
      return boost::dynamic_pointer_cast<const Vector> (d_vVectors.back())->end();
  }

  inline
  Vector::shared_ptr  MultiVector::getVector ( size_t i ) 
  { 
    return d_vVectors[i]; 
  }

  inline
  size_t  MultiVector::getNumberOfSubvectors () 
  { 
    return d_vVectors.size(); 
  }

  inline
  std::string MultiVector::type() const
  {
    return "MultiVector";
  }

  inline
  MultiVector::~MultiVector()
  {
  }

  inline
  MultiVector::vector_iterator  MultiVector::eraseVector ( vector_iterator c )
  {
    return d_vVectors.erase ( c );
  }
  
  inline
  Vector::shared_ptr  MultiVector::create ( Variable::shared_ptr name , AMP_MPI comm )
  {
    Vector::shared_ptr retval ( new MultiVector ( name ) );
    retval->castTo<MultiVector>().d_Comm = comm;
    return retval;
  }

  inline
  AMP_MPI MultiVector::getComm() const
  {
    return d_Comm;
  }
  
  inline
  void MultiVector::dataChanged () 
  { 
    fireDataChange(); 
  }

  inline
  Vector::shared_ptr  MultiVector::create ( const std::string &name , AMP_MPI comm )
  {
    Vector::shared_ptr retval ( new MultiVector ( Variable::shared_ptr ( new MultiVariable ( name ) ) ) );
    retval->castTo<MultiVector>().d_Comm = comm;
    return retval;
  }

  inline
  MultiVector::MultiVector ( Variable::shared_ptr name )
  {
    setVariable ( name );
    d_vGlobalOffsets.push_back ( 0 );
    d_vLocalOffsets.push_back ( 0 );
  }

  inline
  MultiVector::vector_iterator MultiVector::beginVector() 
  { 
    return d_vVectors.begin(); 
  }

  inline
  MultiVector::vector_iterator MultiVector::endVector()   
  { 
    return d_vVectors.end(); 
  }


  inline
  const Vector::shared_ptr  &MultiVector::getVector ( const VectorOperations &rhs , size_t which )  const
  {
    AMP_ASSERT ( which < rhs.castTo<MultiVector>().d_vVectors.size() );
    return rhs.castTo<MultiVector>().d_vVectors[which];
  }
  
  inline
  Vector::shared_ptr  &MultiVector::getVector ( VectorOperations &rhs , size_t which )  const
  {
    AMP_ASSERT ( which < rhs.castTo<MultiVector>().d_vVectors.size() );
    return rhs.castTo<MultiVector>().d_vVectors[which];
  }
  
}
}

