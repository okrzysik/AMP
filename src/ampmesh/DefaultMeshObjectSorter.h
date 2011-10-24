#ifndef included_AMP_DefaultMeshObjectSorter_h
#define included_AMP_DefaultMeshObjectSorter_h

#include "boost/shared_ptr.hpp"
#include "MeshManager.h"

namespace AMP { 
namespace Mesh {


  /** \class DefaultMeshObjectSorter
    * \brief  This class will present the identity permutation for sorting of objects
    * in a mesh.  That is, the mesh object identifiers do not affect the sorting order
    * \tparam ITERATOR  The mesh object iterator this sorter will sort
    */
  template <typename ITERATOR>
  class DefaultMeshObjectSorter : public AMP::LinearAlgebra::ObjectSorter
  {
    protected:
      /** \brief The number of elements in the sorter
        */
      size_t   d_Size;

    public:
      /** \brief  Convenience typedef
        */
      typedef boost::shared_ptr<AMP::LinearAlgebra::ObjectSorter>                  shared_ptr;

      /** \brief Constructor
        * \param[in]  start   The first object to sort
        * \param[in]  end     The last object to sort
        */
      DefaultMeshObjectSorter ( ITERATOR start , ITERATOR end );

      /** \brief Destructor
        */
      virtual ~DefaultMeshObjectSorter ();

      virtual int     operator [] ( int ) const;
      virtual size_t  size () const;
  };

}
}

#include "DefaultMeshObjectSorter.tmpl.h"
#endif
