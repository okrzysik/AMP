#ifndef included_AMP_MultiVariable_h
#define included_AMP_MultiVariable_h

#include "Variable.h"
#include <vector>

namespace AMP {
namespace LinearAlgebra {


/** \brief  A class for combining variables.
 * \details  When physics are brought together, individual variables need
 * to be combined to generate a composition.  For instance, combining
 * temperature and displacement into a single variable.
 *
 * \see MultiVector
 */
class MultiVariable : public Variable
{
public:
    /** \brief Get the first variable in the MultiVariable
     * \return An iterator pointing to the first variable
     */
    inline auto beginVariable() { return d_vVariables.begin(); }

    /** \brief Get end of the MultiVariable array
     * \return An iterator pointing to the end
     */
    inline auto endVariable() { return d_vVariables.end(); }

    /** \brief Get the first variable in the MultiVariable
     * \return An iterator pointing to the first variable
     */
    inline auto beginVariable() const { return d_vVariables.begin(); }

    /** \brief Get end of the MultiVariable array
     * \return An iterator pointing to the end
     */
    inline auto endVariable() const { return d_vVariables.end(); }

    /** \brief If there are multiple matching variables in the list, this
     *  will remove them.  Note that may change the etnry order and will remove any null entries.
     */
    void removeDuplicateVariables();

    /** \brief Given a vector of strings, this will sort the MultiVariable
     * to the given order
     * \param[in] v A list of names by which to sort the MultiVariable
     */
    void sortVariablesByName( const std::vector<std::string> &v );

    /** \brief Constructor
     * \details Because a MultiVariable is a Variable, it must have a name.  This does
     *    not change the names of the variables in the list of vectors.
     * \param[in] name  The name of the MultiVariable
     * \param[in] vars  Optional list of variables in the MultiVariable
     *
     */
    explicit MultiVariable( const std::string &name,
                            const std::vector<std::shared_ptr<Variable>> &vars =
                                std::vector<std::shared_ptr<Variable>>() );

    /** \brief Destructor
     *
     */
    virtual ~MultiVariable();

    /** \brief  Get a particular variable from the list of variables
      * \param  which  the index of the variable sought
      *
      * \details This is an alias for \code
        d_vVariables[which];
        \endcode It is bounds checked in
      * debug builds.
      */
    virtual std::shared_ptr<Variable> getVariable( size_t which );

    /** \brief  Get a particular variable from the list of variables
      * \param  which  the index of the variable sought
      *
      * \details This is an alias for \code
        d_vVariables[which];
        \endcode It is bounds checked in
      * debug builds.
      */
    virtual std::shared_ptr<const Variable> getVariable( size_t which ) const;

    /** \brief Returns the number of variables in the list
      *
      * \details This is an alias for
        \code
        d_vVariables.size();
        \endcode
      */
    virtual size_t numVariables() const;

    /** \brief Add a variable to the end of the variable list
      * \param  newVar  a shared pointer to the new variable
      *
      * \details This is an alias for
        \code
        d_vVariables.push_back ( newVar );
        \endcode
        unless newVar is a MultiVariable.  In order to keep
        heirarchies to a minimum, the members of newVar are added
        instead of newVar itself.
      */
    virtual void add( std::shared_ptr<Variable> newVar );

    /** \brief Set a particular variable in the list
      * \param i    index into the list
      * \param var  a shared pointer to the variable to be placed in the list
      *
      * \details  This is an alias for
        \code
        d_vVariables[i] = var;
        \endcode
        * This is bounds checked in debug builds
      */
    virtual void setVariable( size_t i, std::shared_ptr<Variable> &var );

    // These are adequately documented elsewhere.
    virtual bool operator==( const Variable &rhs ) const override;
    virtual std::shared_ptr<Variable> cloneVariable( const std::string &name ) const override;
    virtual void setUnits( const Units &units ) override;

protected:
    //! List of variables comprising the MultiVariable
    std::vector<std::shared_ptr<Variable>> d_vVariables;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
