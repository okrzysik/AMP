#ifndef  included_AMP_SubsetVariable_H
#define  included_AMP_SubsetVariable_H

#include "VectorIndexer.h"

namespace AMP {
namespace LinearAlgebra {

  /** \class SubsetVariable
    * \brief A variable used to create a SubsetVector
    * \see SubsetVector
    * \see VectorIndexer
    */
  class SubsetVariable : public Variable
  {
    public:
      /** \brief Constructor
        * \param[in]  name  The name of the variable
        */
      SubsetVariable ( const std::string &name );

      /** \brief Return a VectorIndexer that describes the subset
        * \return The VectorIndexer
        */
      virtual VectorIndexer::shared_ptr  getIndexer () = 0;
  };

}
}

#include "SubsetVariable.inline.h"

#endif
