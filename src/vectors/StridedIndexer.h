#ifndef included_AMP_StridedIndexer
#define included_AMP_StridedIndexer

#include "VectorIndexer.h"

namespace AMP {
namespace LinearAlgebra {


  /** \class StridedIndexer
    * \brief A reinterpretation mapping for a vector that selects every \f$i\f$th
    * element in a vector offset by \f$j\f$.
    * \details  A vector of displacements \f$\mathbf{d} = \{x_0 y_0 z_0 x_1 
    * y_1 z_1 \ldots x_n y_n z_n\}^t\f$.  A strided indexer allows a user to
    * extract individual directions such as \f$\mathbf{y} = \{y_0 y_1 \ldots
    * y_n\}^T\f$.
    *
    * Generally, a user will not need to interact with this class.  It can
    * be generated automatically.
    *
    * \see SubsetVector
    * \see VS_Stride
    */
  class  StridedIndexer : public VectorIndexer
  {
    private:
      size_t   d_Offset;
      size_t   d_Stride;

    public:
      /** Constructor
        * \param[in]  offset  The offset from which to stride
        * \param[in]  stride  The number of elements to skip
        */
      StridedIndexer ( size_t  offset ,   size_t  stride );

      virtual bool     isInSub ( size_t ) const;

      virtual size_t   getSubID ( size_t ) const;

      virtual size_t   getSuperID ( size_t ) const;

      virtual size_t   getIncrement ( size_t ) const;

      virtual size_t   getNumLocalElements ( Vector::shared_ptr ) const;
  };

}
}

#include "StridedIndexer.inline.h"

#endif
