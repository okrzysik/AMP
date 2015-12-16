#ifndef included_AMP_DataChangeFirer_h
#define included_AMP_DataChangeFirer_h

#include "utils/shared_ptr.h"
#include <algorithm>
#include <deque>


namespace AMP {
namespace LinearAlgebra {

class DataChangeListener;

/**
  * \class  DataChangeFirer
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
  * The DataChangeFirer provides a mechanism to register and de-register itself with a
  * DataChangeFirer.  It is a list of DataChangeFirers with which it is registered.  It
  * implements a dataChanged() method used by the DataChangeFirer to indicate managed data
  * has changed.
  */

class DataChangeFirer : private std::deque<DataChangeListener *> {
public:
    /**
      * \brief  Notify all listeners of a data change event
      *
      * \details Iterates through the deque of DataChangeListeners and
      * invokes the dataChanged() method.
      */
    virtual void fireDataChange();

    /**
      * \brief  Cosntruct the DataChangeFirer
      *
      * \details When constrcuted, a DataChangeFirer object is
      * empty.  It has no default listeners.
      */
    DataChangeFirer();

    /**
      * \brief  Destroy the DataChangeFirer
      *
      * \details  On destruction, the DataChangeFirer will deregister itself
      * with all DataChangeListeners it has.
      */
    virtual ~DataChangeFirer();

    /**
      * \brief  Register a listener with this DataChangeFirer
      * \param[in]  listener The listener to be alerted to data change events
      * \details Adds a listener to itself.
      */
    virtual void registerListener( DataChangeListener *listener );

    /**
      * \brief  Deregister a listener with this DataChangeFirer
      * \param[in]  listener  The listener to remove from the list
      * \details Removes a listener from itself.
      */
    virtual void deregisterListener( DataChangeListener *listener );
};
}
}


#include "DataChangeFirer.inline.h"

#endif
