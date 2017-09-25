#ifndef included_AMP_MathExpr
#define included_AMP_MathExpr


#include <string>
#include <vector>

#include "utils/tinyexpr/tinyexpr.h"


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
    explicit MathExpr( const std::string &expression,
                       const std::vector<std::string> &variables = std::vector<std::string>() );

    //! Copy constructor
    MathExpr( const MathExpr & );

    //! Assignment operator
    MathExpr &operator=( const MathExpr & );

    //! Move constructor
    MathExpr( MathExpr &&rhs );

    //! Move assignment operator
    MathExpr &operator=( MathExpr &&rhs );

    //! Destructor
    ~MathExpr();

    /**
     * \brief Evaluate
     * \details  Evaluate the expression provided in the constructor with given arguments
     * \param[in] data              List of variable values: { x, y }
     * @return                      Returns the result
     */
    double eval( const std::vector<double> &data = std::vector<double>() );

private:
    void initialize( const std::string &expression, const std::vector<std::string> &variables );

    std::string d_expr;
    te_expr *d_fun;
    std::vector<std::string> d_vars;
    std::vector<te_variable> d_tevar;
    std::vector<double> d_data;
};


} // namespace AMP

#endif
