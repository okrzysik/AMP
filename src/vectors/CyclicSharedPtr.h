#ifndef included_AMP_CyclicSharedPtr
#define included_AMP_CyclicSharedPtr


namespace AMP {
namespace LinearAlgebra {

  /**
    * \class CyclicSharedPtr
    * \brief A class prevents a shared_ptr from being destroyed if it
    * falls out of scope so that if the data comes back into scope, the
    * shared_ptr mechanism can be used.
    * \details
    *  When a vector is created through a TPL clone interface (such as VecClone
    *  or NVClone), the clone created is a shared_ptr.  However, the TPL does
    *  not have a facility for managing a shared_ptr.  When this clone is passed
    *  back to AMP, the shared_ptr mechanism is not set up properly.  By having
    *  the vector hold a shared_ptr to itself, the Vector can be used within AMP.
    *
    *  To use, extend the class which provides createCyclicSharedPtr() and 
    *  destroyCycle().
    */
  template <typename T , typename DELETER>
  class CyclicSharedPtr
  {
    private:
      AMP::shared_ptr<T>  *d_pThis;

    public:
      /** \brief  Construct the CyclicSharedPtr class.  This does not create a cycle.
        */
      CyclicSharedPtr ();

      /** \brief  Destroy the CyclicSharedPtr class.  This requires any created cycle to be previously destroyed.
        */
      virtual ~CyclicSharedPtr ();

      /** \brief  Create a cycle in the shared_ptr for use later.  This should
        *         be called in the TPL method re-implementation of clone.
        */
      void   createCyclicSharedPtr ();

      /** \brief  Destroy the cycle.  This should be called in the TPL method
        *         re-implementation of destroy.
        */
      void   destroyCycle ();
  };

}
}

#include "CyclicSharedPtr.inline.h"
#endif
