#ifndef  included_AMP_DOFMap_h
#define  included_AMP_DOFMap_h

#include "boost/shared_ptr.hpp"
#include "vectors/VectorEntryMap.h"
#include "MeshObject.h"

namespace AMP { 
namespace Mesh {

  /** \class DOFMap
    * \brief A helper class for VectorEntryMap that implements the libMesh interface for
    * mapping mesh objects to entries in a vector
    */
  class DOFMap : public AMP::LinearAlgebra::VectorEntryMap<true>
  {
    public:
      /** \brief  Convenience typedef
        */
      typedef  boost::shared_ptr<DOFMap>     shared_ptr;

      /** \brief Create a DOFMap from a set of parameters.
        * \param[in] rhs  The parameters used to create the DOFMap
        */
      DOFMap ( AMP::LinearAlgebra::VectorEntryMap<true>::Parameters::shared_ptr  rhs );

      /** \brief Get the entry indices of nodal values given an element
        * \param[in]  obj  The element to collect nodal objects for
        * \param[out] ids  The entries in the vector associated with D.O.F.s on the nodes
        * \param[in]  which  Which D.O.F. to get.  If not specified, return all D.O.F.s
        * \details  This will return a vector of pointers into a Vector that are associated with
        * which.
        */
      void getDOFs ( const MeshObject &obj , std::vector <unsigned int> &ids , unsigned int which = static_cast<unsigned int> ( -1 ) ) const;

      /** \brief Get the entry indices of nodal values given a node
        * \param[in]  obj  The node to retrieve ids for
        * \param[out] ids  The indices
        * \param[in]  which  A vector of D.O.F.s needed for the node.  If empty, all indices are returned
        */
      void getDOFs ( const MeshObject &obj , std::vector <unsigned int> &ids , std::vector <unsigned int> which ) const;

      /** \brief  The first D.O.F. on this core
        * \return The first D.O.F. on this core
        */
      size_t  beginDOF ();

      /** \brief  One past the last D.O.F. on this core
        * \return One past the last D.O.F. on this core
        */
      size_t  endDOF ();
  };

}
}

#include "DOFMap.inline.h"

#endif
