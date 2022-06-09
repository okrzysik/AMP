#ifndef included_AMP_MathExpr
#define included_AMP_MathExpr


#include <initializer_list>
#include <memory>
#include <string>
#include <vector>

#include "AMP/utils/tinyexpr/tinyexpr.h"


namespace AMP {


/**
 * \class MathExpr
 * \brief A class used to evaluate math expressions
 * \details  This class provides the ability to evaluate a string
 *     as a math expression.
 *     Note: currently this class is not thread-safe.  Each thread must have a different copy.
 */
class MathExpr
{
public:
    //! Empty constructor
    MathExpr();

    /**
     * \brief Default constructor
     * \details  Construct a MathExpr object to evaluate an expression with input variables
     * \param[in] expression        Expression to evaluate: "sqrt(x^2+y^2)"
     * \param[in] variables         List of variables: { "x", "y" }
     */
    explicit MathExpr( std::string expression,
                       std::vector<std::string> variables = std::vector<std::string>() );

    //! Copy constructor
    MathExpr( const MathExpr & ) = delete;

    //! Assignment operator
    MathExpr &operator=( const MathExpr & ) = delete;

    //! Move constructor
    MathExpr( MathExpr &&rhs );

    //! Move assignment operator
    MathExpr &operator=( MathExpr &&rhs );

    //! Destructor
    ~MathExpr();

    //! Return the expression
    inline const auto &getExpr() const { return d_expr; }

    //! Return the variables
    inline const auto &getVars() const { return d_vars; }

    /**
     * \brief Evaluate
     * \details  Evaluate the expression provided in the constructor with given arguments
     * \param[in] data              List of variable values: { x, y }
     * @return                      Returns the result
     */
    double operator()( const std::initializer_list<double> &data = {} ) const;

    /**
     * \brief Evaluate
     * \details  Evaluate the expression provided in the constructor with given arguments
     * \param[in] data              List of variable values: { x, y }
     * @return                      Returns the result
     */
    double operator()( const double *data ) const;

    //! Check if two expressions are the same
    bool operator==( const MathExpr &rhs ) const;

    //! Check if two expressions are different
    bool operator!=( const MathExpr &rhs ) const { return !operator==( rhs ); }

private:
    void initialize();

    std::string d_expr;
    te_expr *d_fun;
    std::vector<std::string> d_vars;
    std::vector<te_variable> d_tevar;
    mutable std::vector<double> d_data;
};


} // namespace AMP

#endif
