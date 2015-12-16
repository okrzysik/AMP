#ifndef  included_AMP_SubsetVariable_H
#define  included_AMP_SubsetVariable_H

#include "vectors/VectorIndexer.h"
#include "vectors/Variable.h"
#include "discretization/subsetDOFManager.h"


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
      explicit SubsetVariable ( const std::string &name );

      /** \brief Return a DOFManager that describes the subset
        * \return The DOFManager
        * \param[in]  manager  The DOF manager we want to subset
        */
      virtual AMP::Discretization::DOFManager::shared_ptr  getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr manager ) const = 0;
  };

}
}

#include "SubsetVariable.inline.h"

#endif
