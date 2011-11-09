#ifndef included_AMP_DataChangeListener_h
#define included_AMP_DataChangeListener_h

#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <deque>

namespace AMP {
namespace LinearAlgebra {


  class DataChangeFirer;

  /**
    * \class  DataChangeListener
    * \brief  Interface for managing classes that need to know when managed data has changed
    *
    * \details  Some linear algebra packages, such as PETSc, use sophisticated caching
    * schemes to reduce the amount of communication necessary to perform some computations.
    * For instance, unless the values of a vector change, the L1 norm will remain constant.
    * Since computing the L1 norm requires communication for parallel vectors, communication
    * can be reduced by caching this value and invalidating the cache when the data changes.
    * To pass these messages around, a set of classes are provided which allow for the invalidation
    * of the cache.
    * These classes are used when AMP managed data is changed by some linear algebra packages
    * to inform other packages that their cache needs to be invalidated.  These classes
    * implement a callback mechansim in an environment of reference counted pointers.
    *
    * The DataChangeListener provides a mechanism to register and de-register itself with a 
    * DataChangeFirer.  It is a list of DataChangeFirers with which it is registered.  It
    * implements a dataChanged() method used by the DataChangeFirer to indicate managed data
    * has changed.
    */

  class DataChangeListener: private std::deque<DataChangeFirer *>
  {
    private:
      /**
        * \brief  Register this DataChangeListener with a DataChangeFirer
        * \param  firer  the data change firer this object will listen for
        * \details  This method will push firer onto the deque.  This method
        * is only used by the DataChangeFirer
        */
      virtual void        registerWithFirer ( DataChangeFirer *firer );

      /**
        * \brief  Deregister this DataChangeListener with a DataChangeFirer
        * \param  firer  the data change firer this object will no longer listen for
        * \details  This method will erase firer from the deque.  This method
        * is only used by the DataChangeFirer
        */
      virtual void        deregisterFromFirer ( DataChangeFirer *firer );

      friend class DataChangeFirer;

    public:
      /**
        * \brief  Construct the DataChangeListener
        *
        * \details When constructed, a DataChangeListener object is
        * empty.  It has no default firers.
        */
      DataChangeListener ();

      /**
        * \brief  Destroy the DataChangeListener
        *
        * \details On destruction, the DataChangeListener will deregister
        * itself with all DataChangeFirers it is registered with.
        */
      virtual ~DataChangeListener ();

      /**
        * \brief  The method called when a data change event occurs
        *
        * \details  The method called when a data change event occurs
        */
      virtual void        dataChanged () = 0;

  };

}
}

#include "DataChangeListener.inline.h"
#endif
