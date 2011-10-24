#ifndef  included_AMP_StridedVariable_H
#define  included_AMP_StridedVariable_H

#include "StridedIndexer.h"
#include "SubsetVariable.h"

namespace AMP {
namespace LinearAlgebra {

  /** \class StridedVariable
    * \brief An AMP Variable that describes how to stride a vector to create
    * a SubsetVector
    * \see SubsetVector
    * \see StridedIndexer
    */

  class StridedVariable : public SubsetVariable
  {
    private:
      VectorIndexer::shared_ptr   d_Indexer;
      StridedVariable ();

    public:
      /** \brief Constructor
        * \param[in] name  The name of the new variable
        * \param[in] offset  The offset to start striding a vector
        * \param[in] stride  The stride of the vector
        */
      StridedVariable ( const std::string &name , size_t offset , size_t stride );

      virtual VectorIndexer::shared_ptr  getIndexer ();
      virtual  size_t   DOFsPerObject () const;
  };

}
}

#include "StridedVariable.inline.h"

#endif
