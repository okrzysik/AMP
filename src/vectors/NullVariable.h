#ifndef included_AMP_NullVariable_h
#define included_AMP_NullVariable_h


#include "Variable.h"

namespace AMP {
namespace LinearAlgebra {

  /**
    * \brief  A variable that does not require space allocation
    *
    * \details Some variables do not require space to be allocated.  This variable should
    * be used as a placeholder for such variables.
    */

  class NullVariable : public Variable
  {
    public:
      /**
        *  \brief Constructor
        *  \param name the name of the NullVariable
        *
        *  Even though there is no space allocated, this variable still needs a name
        */
      NullVariable ( const std::string &name );

      /**
        * \brief  returns -1
        *
        *  How does a function with an unsigned return give -1?  
        \code
        return static_cast<size_t>(-1);
        \endcode
        *  that's how.
        */
      size_t  variableID () const;

      /**
        * \brief  return number of DOFs per object, which is 0 for NullVariables
        *
        * \details  Even though this function will soon be deprecated, it still
        * returns 0 in this case.
        */
      size_t  DOFsPerObject () const;

      /**
        * \brief  return a clone of this NullVariable
        * \param  name  the name of the new variable
        *
        *  Returns a copy of this variable.
        */
      Variable::shared_ptr  cloneVariable ( const std::string &name ) const;
  };

}
}

#endif
