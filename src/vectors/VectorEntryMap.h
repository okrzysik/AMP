#ifndef included_AMP_VectorEntryMap_h
#define included_AMP_VectorEntryMap_h

#include <map>
#include <vector>

#include "utils/AMP_MPI.h"
#include "utils/ParameterBase.h"
#include "vectors/CommunicationList.h"

#include "ObjectSorter.h"
#include "Variable.h"

namespace AMP {
namespace LinearAlgebra {

  /** \brief Parameters used to instantiate a VectorEntryMap
    */
  class VectorEntryMapParameters 
  {
    public:
      /** \brief Convenince typedef 
        */
      typedef AMP::shared_ptr<VectorEntryMapParameters>    shared_ptr;

      /** \brief A sorting of the entities of the vector
        */
      ObjectSorter::shared_ptr          d_SortedObjects;

      /** \brief A variable describing the vector
        */
      Variable::shared_ptr              d_Variable;

      /** \brief The communicator to construct the VectorEntryMap over
        */
      AMP_MPI                           d_Comm;

      /** \brief The communication list of the vector
        */
      CommunicationList::shared_ptr     d_CommList;
  };

  /** \class VectorEntryMap
    * \tparam AFFINE_MAP  Whether the mapping is an affine map characterized by \f$\alpha\f$ 
    * and \f$\beta\f$ as  \f$\mathrm{entry}_\mathrm{id} = \alpha \mathrm{object}_\mathrm{id} + 
    * \beta \f$.
    * \brief A mapping from an object ID to a global ID
    * \details  A Vector presents a contiguous block of data indexed from 0 on each core.
    * This class will associate a unique global ID across an MPI Communicator to each
    * object provided in the ObjectSorter.  The global IDs are packed, entries \f$[0\ldots n_1-1]\f$ 
    * are on the first core, \f$[n_1\ldots n_2-1]\f$ are on the second core, etc.
    */
  template <bool AFFINE_MAP>
  class VectorEntryMap : public Castable
  {
    protected:
      /** \brief The CommunicationList for this mapping
        */
      CommunicationList::shared_ptr     d_CommList;

      /** \brief The sort of the objects used for this mapping
        */
      ObjectSorter::shared_ptr          d_SortedObjects;

      /** \brief The variable used to describe this mapping
        */
      Variable::shared_ptr              d_Variable;

      /** \brief The non-affine map if it exists
        */
      std::map<size_t,size_t>           d_NonAffineMap;

      /** \brief The communicator for this mapping
        */
      AMP_MPI                           d_Comm;

      /** \brief  The number of local rows in the vector
        */
      size_t                            d_NumRows;

      /** \brief  The first d.o.f. in the vector
        */
      size_t                            d_FirstRow;

      /** \brief  The total number of entries in the vector
        */
      size_t                            d_TotalSize;

      /** \brief  Compute the first DOF on this core
        * \return  The first DOF on this core
        */
      size_t                            ComputeFirstDOF ();

    public:
      /** \brief  Convenience typedef
        */
      typedef  VectorEntryMapParameters                   Parameters;

      /** \brief  Convenience typedef
        */
      typedef  AMP::shared_ptr<VectorEntryMap>          shared_ptr;

      /** \brief Constructor
        * \param[in] ptr  The parameters used to create this VectorEntryMap
        */
      VectorEntryMap ( Parameters::shared_ptr ptr );

      /** \brief  Destructor
        */
      virtual ~VectorEntryMap();


      /** \brief  Return the global id in the vector given an object and associated dof
        * \param[in]  obj_id   The id of the object
        * \param[in]  dof   The d.o.f. of the object with respect to the object
        * \details  For instance, given node 7 and a position vector, this can be used to find
        * the element in the vector for the \f$y\f$-coordinate if node 7
        \code
        size_t  y_coord = disp.getGlobalID ( 7 , 1 );
        \endcode
        */
      size_t getGlobalID ( size_t obj_id , size_t dof ) const;

      /** \brief  Given an element in the vector, return the object id and dof
        * \param[in] global  The vector global id
        * \param[out] obj_id  The object associated with that global id
        * \param[out] dof  The d.o.f. associated with that global id
        */
      void   getLocalID ( size_t global , size_t &obj_id , size_t &dof ) const;

      /** \brief  Build a non-affine map from an iterator
        * \tparam DOF_FIELD  A function that returns the number of D.O.F. given an iterator
        * \tparam ITERATOR  The iterator type to search over
        * \param[in] begin  The start of the list to compute the non-affine map
        * \param[in] end    The end of the list to compute the non-affine map
        * \details  This will create a non affine map, given an integer field over the set
        * provided by the iterators.  For instance, to create 3 d.o.f. for each of a NodeIterator:
        \code
        template <size_t  I , typename ITERATOR>
        size_t   constDofField ( const ITERATOR )
        {
          return I;
        }
         .
         .
         .
       NodeSetIterator  firstNode = NodeSet.begin();
       NodeSetIterator  endNode = NodeSet.end();
       VectorEntryMap<false>  mapping (...);
       mapping.computeNonAffineMap<constDofField<3,NodeSetIterator>,NodeSetIterator> ( firstNode, endNode);
       \endcode
       */
      template <typename DOF_FIELD , typename ITERATOR>
      void  computeNonAffineMap ( ITERATOR begin , ITERATOR end );

      /** \brief  Get the communication list for this mapping
        * \return  The CommunicationList
        */
      CommunicationList::shared_ptr   &getCommunicationList ();

      /** \brief  Set the communication list for this mapping
        * \param[in] in  The CommunicationList
        */
      void  setCommunicationList ( CommunicationList::shared_ptr in );

      /** \brief Return the number of DOFs per object
        * \return  d_Variable->DOFsPerObject()
        */
      virtual  unsigned int   numDOFsPerObject ();

      /** \brief Return the number of local elements
        * \return  the number of entries on this core
        */
      virtual  unsigned int   numLocalElements ();

      /** \brief  Return the first element on this core
        * \return  the global id of the first entry on this core
        */
      virtual  unsigned int   firstElement ();

      /** \brief  Return one more than the last element on this core
        * \return  last global id + 1
        */
      virtual  unsigned int   endElement ();
      
      /** \brief  Return total number of global elements
        * \return  the number of indexed elements
        */
      virtual  unsigned int   numGlobalElements ();
  };

}
}

#include "VectorEntryMap.tmpl.h"

#endif
