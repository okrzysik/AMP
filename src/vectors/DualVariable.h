#ifndef included_AMP_DualVariable_h
#define included_AMP_DualVariable_h

#include "Variable.h"

namespace AMP {
namespace LinearAlgebra {

  /**
    * \class DualVariable
    * \brief A composition of two variables into one
    * \details  For an operator \f$\mathcal{F} (\mathbf{x}, T)\f$,
    * there are two variables: a displacement \f$\mathbf{x}\f$, and a
    * temperature \f$T\f$.  These two variables can be combined into one
    * for the purpose of describing the input of the operator 
    * \f$\mathcal{F}\f$.
    *
    * \see DualVector
    */
  class DualVariable : public Variable
  {
    protected:

      /**
        * \brief The first variable in the composition
        */
      Variable::shared_ptr   d_Var1;

      /**
        * \brief The second variable in the composition
        */
      Variable::shared_ptr   d_Var2;

    public:
      /**
        * \brief Construct a DualVariable from two variables and a name
        * \param[in] one The first variable
        * \param[in] two The second variable
        * \param[in] name The name of the dual variable
        */
      DualVariable ( Variable::shared_ptr one , Variable::shared_ptr two , std::string name = "noname" );

      /**
        * \brief Destroy the dual variable
        */
      virtual ~DualVariable ();

      /**
        * \brief Get the first variable
        * \return  A shared pointer to the first variable
        */
      Variable::shared_ptr  first ();

      /**
        * \brief Get the second variable
        * \return  A shared pointer to the second variable
        */
      Variable::shared_ptr  second ();

      /**
        * \brief Return the hash value that a mesh uses to cache
        * vectors, matrices, and DOFmaps
        * \warning  Always throws an exception
        * \details  This method throws an execption.  Rather,
        * the variableID of the constituent variables should be used.
        */
      virtual  size_t  variableID () const;

      /**
        * \brief Compare this variable to another
        * \param[in] rhs The variable to comaper with
        * \return true if variables are equal
        * \details Returns true iff rhs is a DualVariable and d_Var1 = rhs.d_Var1 and d_Var2 = rhs.d_Var2
        */
      virtual bool   operator == ( const Variable &rhs ) const;

      /**
        * \brief Create a clone of this variable
        * \param[in] name The name given to the cloned variable
        * \details This will copy the shared pointers of the constituent
        * variables to the clone.
        * \return A shared pointer to a clone of this variable
        */
      virtual Variable::shared_ptr  cloneVariable ( const std::string &name ) const;

      /**
        * \brief The number of DOFs per object for this variable
        * \details This throws an exception for DualVariables.
        * \warning Always throws an exception
        */
			virtual  size_t   DOFsPerObject () const;
  };

}
}

#include "DualVariable.inline.h"

#endif
