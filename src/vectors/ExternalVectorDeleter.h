#ifndef included_AMP_ExternalVectorDeleter
#define included_AMP_ExternalVectorDeleter

namespace AMP {
namespace LinearAlgebra {

  class Vector;

  /** \brief Class used to prevent the shared pointer from deleting a Vector
    * unexpectedly
    * \details  In certain situations, it is necessary to have a shared pointer
    * of a Vector without the shared pointer having control over memory management.
    * This class provides this functionality
    */
  class ExternalVectorDeleter
  {
    public:
      /** \brief Empty delete method
        * \param[in] v  Vector not to delete
        */
      void operator()(Vector *v);
  };

}
}

#include "ExternalVectorDeleter.inline.h"
#endif
