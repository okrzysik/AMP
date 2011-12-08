#ifndef included_AMP_ManagedDataMap
#define included_AMP_ManagedDataMap

#include <vector>
#include <stdexcept>

namespace AMP {
namespace LinearAlgebra {

/** \brief For libraries that do not assume "owned" vector elements are numbered
  * contiguously, this class provides a non-contiguous mapping
  * \details  This creates
  */
class ManagedDataMap 
{
public:
    typedef std::vector<size_t>                  mapping;
    typedef std::vector<size_t>::iterator        iterator;
    typedef std::vector<size_t>::const_iterator  const_iterator;
    typedef boost::shared_ptr<mapping>           mapping_ptr;


    //! Copy constructor
    ManagedDataMap ( const ManagedDataMap & );

    /** \brief Constructor
      * \param[in] local_size  The number of vector entries on this core
      * \param[in] global_size  The number of vector entries on all cores
      */
    ManagedDataMap ( size_t local_size , size_t global_size );

    /** \brief Map a local row id to a global row id
      * \param[in] local_id  The local id to be associated
      * \param[in] global_id  The global id to be associated
      * \details  This will associate a local id to a global id
      */
    void  addMapping ( size_t local_id , size_t global_id );

    /** \brief  Return the global_id of the first row of the map.
      * \return  The global_id
      */
    size_t   firstRow ();

    /** \brief Return the iterator pointing to the first entry.
      * \return An iterator to the first entry in the map
      */
    iterator   begin();

    /** \brief Return the end iterator
      * \return The end iterator
      */
    iterator   end();

    /** \brief Return the iterator pointing to the first entry.
      * \return An iterator to the first entry in the map
      */
    const_iterator  begin() const; 
    /** \brief Return the end iterator
      * \return The end iterator
      */
    const_iterator  end() const;

    /** \brief Get the number of entries in the map on this core
      * \return  The number of entries
      */
    size_t   getLocalSize () const;

    /** \brief Get the number of entries in the map on all cores
      * \return  The number of entries
      */
    size_t   getGlobalSize () const;

protected:
    //!  Blank constructor
    ManagedDataMap () {}

    //!  The global id of the first row in the local vector
    size_t                      d_FirstRow;

    //!  The mapping
    mapping_ptr                 d_vMapping;

    //!  Number of local entries
    size_t                      d_uiLocalSize;
      
    //!  Number of global entries
    size_t                      d_uiGlobalSize;

    bool                        d_FirstRowFilled;

};


}
}

#include "ManagedDataMap.inline.h"
#endif
