#ifndef included_AMP_SubsetVector
#define included_AMP_SubsetVector

#include "AMP/vectors/Vector.h"
#include "AMP/vectors/operations/VectorOperationsDefault.h"
#include <vector>

namespace AMP {
namespace LinearAlgebra {

/** \class SubsetVector
  * \brief This vector is a subset of an AMP Vector
  * \details
    Given an AMP Vector, this class will create a view of a subset of the
    vector.  For instance, if \f$\mathbf{a} = \{ a_0 a_1 a_2 \ldots a_n\}\f$,
    and \f$S\f$ is a set of non-negative integers \f$( s_0 s_1 \ldots s_m )\f$
    such that \f$0 \le s_i \le n\f$, then the SubsetVector will be the
    vector \f$\mathbf{a}_S = \{ a_{s_0} a_{s_1} a_{s_2} \ldots a_{s_m}\} \f$.

    This class provides a factory method called view:
    \code
      AMP::LinearAlgebra::Vector::shared_ptr  vec1;
      AMP::LinearAlgebra::Vector::shared_pt   vec2 = AMP::SubsetVector( vec1 , subsetVar );
      AMP::LinearAlgebra::Vector::shared_ptr  vec3 = vec2->clone( "subset2" );
    \code

    Since this is a view, any change to vec2 will be reflected on vec1 and
    vice versa.  vec2 is a sparse vector, with a mapping of new index to old.
    vec3, on the other hand, is a dense vector without an index.  If a lot
    of computation is necessary on the sparse vector, the data can be copied
    in and out:

    \code
      // Subset the vector to make a sparse vector.
      vec2 = AMP::SubsetVector( vec1 , subsetVar )

      // Copy the sparse vector data to a dense vector
      vec3.copyVector( vec2 );

      // Perform whatever math
      performComputation( vec3 );

      // Copy data back to vec2, and, consequently, vec1
      vec2.copyVector( vec3 );
    \endcode
  */
class SubsetVector : public Vector
{

public:
    static Vector::shared_ptr view( Vector::shared_ptr, Variable::shared_ptr );
    static Vector::const_shared_ptr view( Vector::const_shared_ptr, Variable::shared_ptr );
    std::string type() const override;
    using Vector::cloneVector;
    Vector::shared_ptr cloneVector( Variable::shared_ptr ) const override;
    void swapVectors( Vector &rhs ) override;
    void aliasVector( Vector &rhs ) override;
    void assemble() override {}
    
private:
    SubsetVector() {}
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
