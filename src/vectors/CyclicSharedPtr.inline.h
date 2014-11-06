
#include "utils/Utilities.h"

namespace AMP {
namespace LinearAlgebra {

  template <typename T , typename DELETER>
  CyclicSharedPtr<T,DELETER>::CyclicSharedPtr ()
  {
    d_pThis = 0;
  }

  template <typename T , typename DELETER>
  CyclicSharedPtr<T,DELETER>::~CyclicSharedPtr ()
  {
    AMP_ASSERT ( d_pThis == 0 );
  }

  template <typename T , typename DELETER>
  void CyclicSharedPtr<T,DELETER>::createCyclicSharedPtr ()
  {
    d_pThis = new AMP::shared_ptr<T> ( dynamic_cast<T *> ( this ) , DELETER () );
  }

  template <typename T , typename DELETER>
  void CyclicSharedPtr<T,DELETER>::destroyCycle ()
  {
    delete d_pThis;
    d_pThis = 0;
  }

}
}
