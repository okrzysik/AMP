#ifndef included_AMP_Castable
#define included_AMP_Castable

#include "utils/Utilities.h"

namespace AMP {

/** \class Castable
  * \brief A convenience class for determining type information
  */
class Castable {
public:
    /** \brief Destructor
      */
    virtual ~Castable();

    /** \brief Return a reference to this object of type T
      * \tparam T  The type to cast to
      * \return  A reference to this casted appropriately
      * \details  This method throws an exception if the
      * cast fails
      */
    template <typename T>
    T &castTo();


    /** \brief Return a reference to this object of type T
      * \tparam T  The type to cast to
      * \return  A reference to this casted appropriately
      * \details  This method throws an exception if the
      * cast fails
      */
    template <typename T>
    const T &castTo() const;

    /** \brief Return true if this is of type T
      * \tparam T type to compare against
      * \return  True if this is of type T, false otherwise.
      */
    template <typename T>
    bool isA() const;
};
}

#include "Castable.tmpl.h"

#endif
