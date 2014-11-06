#ifndef included_AMP_VectorIndexer
#define included_AMP_VectorIndexer

#include "utils/shared_ptr.h"
#include "utils/Castable.h"

namespace AMP {
namespace LinearAlgebra {

  class Vector;

  /**
    * \class Indexer
    * \brief  A remapping class used by SubsetVector to translate
    * global indices in a subset vector to the underlying vector
    * \details  A SubsetVector will make a subset of a regular vector
    * look like an AMP vector.  For instance, given a displacement
    * vector, a SubsetVector can be used to extract individual 
    * coordinates.  This class is used by the SubsetVector to translate
    * between global ids.
    * \see SubsetVector
    */
  class VectorIndexer : public Castable
  {
    public:
      typedef AMP::shared_ptr<VectorIndexer>   shared_ptr;

      /** \brief Destructor */
      virtual ~VectorIndexer ();

      /** \brief Given a global id in the underlying vector, this returns true
        * if the global id is in the SubsetVector.
        * \param[in] underlyingID  The global ID of the underlying vector
        * \return True if the underlying ID is in the subset.  False otherwise.
        */
      virtual bool  isInSub ( size_t underlyingID ) const = 0;

      /** \brief Given a global id in the underlying vector, return the subset global id.
        * \param[in] underlyingID  The global ID of the underlying vector
        * \return The global ID in the SubsetVector
        */
      virtual size_t  getSubID ( size_t underlyingID ) const = 0;

      /** \brief Given a global id in the SubsetVector, return the corresponding global id
        * in the underlying vector
        * \param[in] subsetID  The global ID of the SubsetVector
        * \return The global ID in underlying vector
        */
      virtual size_t  getSuperID ( size_t subsetID ) const = 0;

      /** \brief Given a global id in the underlying vector, return the distance to the
        * next id in the SubsetVector.
        * \param[in] underlyingID The global ID of the underlying vector
        * \return getSubID ( underlyingID + 1 ) - getSubID ( underlyingID )
        */
      virtual size_t  getIncrement ( size_t underlyingID ) const = 0;

      /** \brief Given a particular underlying vector, compute the number of local
        * elements that are in this mapping
        * \param[in]  v A vector
        * \return The number of local objects in the vector also in the subset.
        */
      virtual size_t  getNumLocalElements ( AMP::shared_ptr<Vector> v ) const = 0;
  };

}
}

#include "VectorIndexer.inline.h"

#endif
