#ifndef included_AMP_Variable_h
#define included_AMP_Variable_h

#include "utils/Utilities.h"
#include <utils/Castable.h>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace AMP {
namespace LinearAlgebra {

  /**
    * \class Variable
    * \brief A description of the memory allocation required for a vector.
    *
    * \details  An operator \f$D: X\rightarrow Y\f$ transforms a vector in
    * on space into a vector in another.  The operators in AMP work similarly.
    * The vectors in AMP mesh are encapsulated in the Vector class and its
    * derivatives.  For most applications, these vectors are discretizations,
    * a representative of a vector in a function space sampled at points.
    * The vectors can be thought of as a composition of two concepts:  the
    * field that is the range of the function and the discretization.  A
    * Variable encapsulates this field and how it is discretized.  In FEA
    * applications, for example, the field discretization is a mesh.  Given
    * a variable and an explicit discretization, as in a hexahedral mesh of
    * a unit cube, the memory requirements and memory indexing required for
    * simulation can be computed.
    *
    * For instance, if \f$D\f$ requires three degrees of freedom per node,
    * then the variable class encapsulates this information.  When the
    * variable is combined with a particular mesh, a vector of the appropriate
    * size can be created.
    */

  class Variable : public Castable , public boost::enable_shared_from_this<Variable>
  {
    protected:
      /** \brief  A name given to the variable
        *
        * \details Variables have names for easy identification.  For instance,
        * some variables are called displacement, concentration, search direction,
        * etc.
        */
      std::string  d_VariableName;

      /** \brief  The units this variable is measured in
        */
      std::string  d_Units;

    public:
      /** \brief  Set the units of this variable
        */
      virtual void  setUnits ( const std::string &t );

      /** \brief  Get the units of this variable
        */
      virtual const std::string &getUnits () const;

      /** \brief  Shared pointer name
        *
        */
      typedef boost::shared_ptr<Variable>    shared_ptr;

      /** \brief  Construct a variable with a name
        * \param  name  The name of the variable
        *
        * \details At the very least, a variable must have a name.  Since this class
        * is virtual, this constructor is used by derived classes.
        */
      Variable ( const std::string &name );

      /** \brief  Destructor
        *
        */
      virtual ~Variable ();

      /** \brief  A function that returns the name of a variable
        *
        * \details This gives access to the name
        */
      virtual  const std::string  &getName() const;

      /** \brief  Internal method used by meshes to speed memory allocation
        *
        * \details The variableID does not refer to the identifier of the particular
        *  variable but an identifier of the type of memory allocation requested.
        */
      virtual  size_t   variableID () const;

      /** \brief  Compares two variables for equality.
        * \param  rhs  Variable to compare 
        *
        * \details This operation compares not only names but memory allocation pattern.
        * A "temperature" stored for each node is different from a "temperature"
        * stored for each cell
        */
      virtual  bool     operator == ( const Variable & rhs ) const;

      /** \brief  Inverse of ==
        * \param  rhs  Variable to compare
        *
        * \details This function performs an equality check and negates it.  Hence, it
        *  is not virtual
        */
      bool operator != ( const Variable & rhs ) const;

      /** \brief  Set the name of a variable
        * \param  NewName  the new name of the variable
        *
        * \details Due to the hierarchical nature of Vectors and Variables, changing
        *  the name of a variable can be a Very Dangerous Thing.  Before invoking
        *  this function, make sure this is what you really want to do.
        */
      virtual void setName ( const std::string &NewName );

      /** \brief  Create a variable of the same type with a new name
        * \param  name  The name of the new variable
        *
        * \details This function will create a "deep" copy of this variable.
        */
      virtual Variable::shared_ptr   cloneVariable ( const std::string &name ) const;

      /** \brief  Create a variable of the same type with the name "noname"
        *
        * \details This is an alias for Vector::cloneVariable ( "noname" )
        */
      Variable::shared_ptr  cloneVariable () const;

      /** \brief  Return the number of DOFs required per mesh object
        * 
        * \details  This is a convenience function for allocating space.  This
        * function is slated to be deprecated.  DOFsPerObject may not be a constant.
        * Rather, it may be a field on the domain.
        */
      virtual  size_t   DOFsPerObject () const;



      virtual  Variable::shared_ptr  getVariable ( size_t which );
  };

}
}

#include "Variable.inline.h"
#endif
