#ifndef included_AMP_ObjectSorter_h
#define included_AMP_ObjectSorter_h

#include "utils/shared_ptr.h"

namespace AMP {
namespace LinearAlgebra {

/** \class ObjectSorterParameters
  \brief A list of parameters used for creating an ObjectSorter
  \see ObjectSorter
  */
class ObjectSorterParameters : public ParameterBase
{
public:
    /** Convenience typedef
      */
    typedef AMP::shared_ptr<ObjectSorterParameters> shared_ptr;

    /** Constructor
      */
    ObjectSorterParameters();

    /** The id of the first object on a core
      */
    size_t d_FirstObject;
};

/** \class ObjectSorter
  * \brief An abstract class for mapping ids to and from particular orderings
  */
class ObjectSorter
{
private:
    ObjectSorter( const ObjectSorter &rhs );

protected:
    /** Default constructor
      */
    ObjectSorter();

    /** The id of the first object on a core
      */
    size_t d_FirstObject;

public:
    /** Convenience typedef
      */
    typedef AMP::shared_ptr<ObjectSorter> shared_ptr;

    /** Convenience typedef
      */
    typedef ObjectSorterParameters Parameters;

    /** Tag for sort type
      */
    enum class SortType { DefaultNodeOrder };

    /** Create an ObjectSorter from a set of parameters
      */
    ObjectSorter( Parameters::shared_ptr params );

    /** Destructor
      */
    virtual ~ObjectSorter();

    /** \brief Return the first object on this core.
      \return The id of the first object on this core.
      */
    size_t getFirstObject();

    /** \brief Map operation
      * \param[in] id The id in the preimage of the map
      * \return The image of the id under the map
      */
    virtual int operator[]( int id ) const = 0;

    /** \brief  Returns the number of objects in the image on this core
      * \return  The size of the image of the map
      */
    virtual size_t size() const = 0;
};
}
}

#include "ObjectSorter.inline.h"

#endif
