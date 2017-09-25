
#include <algorithm>
#include <math.h>
#include <stdexcept>

#include "utils/Utilities.h"

#include "MultiVector.h"

namespace AMP {
namespace LinearAlgebra {


inline Vector::shared_ptr MultiVector::getVector( size_t i ) { return d_vVectors[i]; }

inline Vector::const_shared_ptr MultiVector::getVector( size_t i ) const { return d_vVectors[i]; }

inline size_t MultiVector::getNumberOfSubvectors() const { return d_vVectors.size(); }

inline std::string MultiVector::type() const { return "MultiVector"; }

inline MultiVector::~MultiVector() {}

inline AMP_MPI MultiVector::getComm() const { return d_Comm; }

inline void MultiVector::dataChanged() { fireDataChange(); }

inline MultiVector::vector_iterator MultiVector::beginVector() { return d_vVectors.begin(); }

inline MultiVector::vector_iterator MultiVector::endVector() { return d_vVectors.end(); }

inline const Vector::shared_ptr &MultiVector::getVector( const VectorOperations &rhs,
                                                         size_t which ) const
{
    auto x = dynamic_cast<const MultiVector *>( &rhs );
    AMP_ASSERT( x != nullptr );
    AMP_ASSERT( which < x->d_vVectors.size() );
    return x->d_vVectors[which];
}

inline Vector::shared_ptr &MultiVector::getVector( VectorOperations &rhs, size_t which ) const
{
    auto x = dynamic_cast<MultiVector *>( &rhs );
    AMP_ASSERT( x != nullptr );
    AMP_ASSERT( which < x->d_vVectors.size() );
    return x->d_vVectors[which];
}


} // namespace LinearAlgebra
} // namespace AMP
