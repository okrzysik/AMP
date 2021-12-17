#ifndef included_AMP_Variable_h
#define included_AMP_Variable_h

#include "AMP/utils/Units.h"

#include <memory>


namespace AMP::LinearAlgebra {


class VectorSelector;


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
     * \details  At the very least, a variable must have a name.
     *    Since this class is virtual, this constructor is used by derived classes.
     * \param  name  The name of the variable
     */
    explicit Variable( const std::string &name );


    //!  Destructor
    virtual ~Variable();


    //!  Set the units of this variable
    virtual void setUnits( const Units &u );


    //!  Get the units of this variable
    virtual const Units &getUnits() const;


    /** \brief  A function that returns the name of a variable
     * \details  This gives access to the name
     */
    virtual const std::string &getName() const;


    /** \brief  Compares two variables for equality.
     * \details This operation compares the names.
     * \param  rhs  Variable to compare
     */
    virtual bool operator==( const Variable &rhs ) const;


    /** \brief  Inverse of ==
     * \details This function performs an equality check and negates it.
     * \param  rhs  Variable to compare
     */
    bool operator!=( const Variable &rhs ) const;


    /** \brief  Create a variable of the same type with a new name
     * \details This function will create a "deep" copy of this variable.
     * \param  name  The name of the new variable
     */
    virtual std::shared_ptr<Variable> cloneVariable( const std::string &name ) const;


    /** \brief  Create a VectorSelector
     * \details This function will create a VectorSelector that is able to subset
     *    a vector for the variable.  This is used by Vector::subsetVectorForVariable
     *    which then calls Vector::selectInto.
     */
    virtual std::shared_ptr<VectorSelector> createVectorSelector() const;


protected:
    /** \brief  A name given to the variable
     * \details Variables have names for easy identification.
     *     For instance, some variables are called displacement, concentration,
     *     search direction, etc.
     */
    std::string d_VariableName;

    /** \brief  The units this variable is measured in
     */
    AMP::Units d_Units;
};

} // namespace AMP::LinearAlgebra

#endif
