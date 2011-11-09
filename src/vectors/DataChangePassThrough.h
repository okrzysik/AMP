#ifndef included_AMP_DataChangePassThrough_h
#define included_AMP_DataChangePassThrough_h

#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <deque>

#include "DataChangeFirer.h"

namespace AMP {
namespace LinearAlgebra {

  /**
    * \class  DataChangePassThrough
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
    * Some Vector types are "wrappers" or "views" of other vectors.  For instance, a MultiVector
    * makes several vectors look like a single vector.  In order for a DataChangeFirer to
    * propogate information to all DataChangeListeners that need to hear about the data change,
    * this class implements both interfaces:  the dataChanged method invokes the fireDataChange.
    *
    * For example, given a PETSc wrapper around a MultiVector of EpetraVectors, the PETSc wrapper
    * needs to invalidate its cache whenever an Epetra method changes one of the vectors of the
    * MultiVector.  In this case, the MultiVector is a DataChangePassThrough.
    *
    * This class will pass along a dataChanged() callback to all registered listeners.
    */

  class DataChangePassThrough : public DataChangeFirer , public DataChangeListener
  {
    public:
      /** \brief Construct the object
        * \details  This is an empty function
        */
      DataChangePassThrough ();

      /** \brief Destroy the object
        * \details  This is an empty function since the base classes handle deregistration
        */
      virtual ~DataChangePassThrough ();

      /** \brief Invoke the dataChange method on all registered listeners.
        * \details  This simply calls the fireDataChange() method of the DataChangeFirer class
        */
      virtual void        dataChanged ();
  };

}
}

#include "DataChangePassThrough.inline.h"

#endif
