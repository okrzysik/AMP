#ifndef included_AMP_Variable_h
#define included_AMP_Variable_h

#include "AMP/utils/Units.h"

#include <memory>


namespace AMP::LinearAlgebra {


/**
 * \class Variable
 * \brief A description of the data in the vector
 *
 * \details  This class stores information about the vector such as the name and units.
 *   It is used in subsetting operations.
 */
class Variable
{
public:
    /** \brief  Construct a variable with a name
     * \param  name  The name of the variable
     *
     * \details At the very least, a variable must have a name.  Since this class
     * is virtual, this constructor is used by derived classes.
     */
    explicit Variable( const std::string &name );


    //!  Destructor
    virtual ~Variable();


    //!  Set the units of this variable
    virtual void setUnits( const Units &u );


    //!  Get the units of this variable
    virtual const Units &getUnits() const;


    /** \brief  A function that returns the name of a variable
     *
     * \details This gives access to the name
     */
    virtual const std::string &getName() const;


    /** \brief  Compares two variables for equality.
     * \param  rhs  Variable to compare
     *
     * \details This operation compares the names.
     * A "temperature" stored for each node is different from a "temperature"
     * stored for each cell
     */
    virtual bool operator==( const Variable &rhs ) const;


    /** \brief  Inverse of ==
     * \param  rhs  Variable to compare
     *
     * \details This function performs an equality check and negates it.  Hence, it
     *  is not virtual
     */
    bool operator!=( const Variable &rhs ) const;


    /** \brief  Create a variable of the same type with a new name
     * \param  name  The name of the new variable
     *
     * \details This function will create a "deep" copy of this variable.
     */
    virtual std::shared_ptr<Variable> cloneVariable( const std::string &name ) const;

protected:
    /** \brief  A name given to the variable
     *
     * \details Variables have names for easy identification.  For instance,
     * some variables are called displacement, concentration, search direction,
     * etc.
     */
    std::string d_VariableName;

    /** \brief  The units this variable is measured in
     */
    AMP::Units d_Units;
};

} // namespace AMP::LinearAlgebra

#endif
