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
      typedef std::vector<int>                  mapping;
      typedef std::vector<int>::iterator        iterator;
      typedef std::vector<int>::const_iterator  const_iterator;
      typedef boost::shared_ptr<mapping>        mapping_ptr;

      bool                         d_FirstRowFilled;

    protected:
      /** \brief The global id of the first row in the local vector
        */
      int                          d_FirstRow;

      /** \brief  The mapping
        */
      mapping_ptr                  d_vMapping;

      /** \brief  Number of local entries
        */
      size_t                       d_uiLocalSize;
      
      /** \brief  Number of global entries
        */
      size_t                       d_uiGlobalSize;

      /** \brief  Blank constructor
        */
      ManagedDataMap () {}

    public:

      /** \brief Copy constructor
        */
      ManagedDataMap ( const ManagedDataMap & );

      /** \brief Constructor
        * \param[in] local_size  The number of vector entries on this core
        * \param[in] global_size  The number of vector entries on all cores
        */
      ManagedDataMap ( int local_size , int global_size );

      /** \brief Map a local row id to a global row id
        * \param[in] local_id  The local id to be associated
        * \param[in] global_id  The global id to be associated
        * \details  This will associate a local id to a global id
        */
      void  addMapping ( int local_id , int global_id );

      /** \brief  Return the global_id of the first row of the map.
        * \return  The global_id
        */
      int   firstRow ();

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
  };

}
}

#include "ManagedDataMap.inline.h"
#endif
