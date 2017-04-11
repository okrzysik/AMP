#include "utils/MathExpr.h"
#include "utils/Utilities.h"
#include "utils/tinyexpr/tinyexpr.h"


namespace AMP {


/************************************************************
*  Constructors/Destructor                                  *
************************************************************/
MathExpr::MathExpr( ):
    d_fun(nullptr)
{
}
MathExpr::MathExpr( const std::string& expression, const std::vector<std::string>& variables )
{
    initialize( expression, variables );
}
MathExpr::MathExpr( const MathExpr &rhs ):
    MathExpr( rhs.d_expr, rhs.d_vars )
{
}
MathExpr& MathExpr::operator=( const MathExpr &rhs )
{
    if ( this == &rhs )
        return *this;
    initialize( rhs.d_expr, rhs.d_vars );
    return *this;
}
MathExpr::MathExpr( MathExpr &&rhs )
{
    std::swap( d_expr,  rhs.d_expr );
    std::swap( d_fun,   rhs.d_fun );
    std::swap( d_vars,  rhs.d_vars );
    std::swap( d_tevar, rhs.d_tevar );
    std::swap( d_data,  rhs.d_data );
}
MathExpr& MathExpr::operator=( MathExpr &&rhs )
{
    if ( this == &rhs )
        return *this;
    std::swap( d_expr,  rhs.d_expr );
    std::swap( d_fun,   rhs.d_fun );
    std::swap( d_vars,  rhs.d_vars );
    std::swap( d_tevar, rhs.d_tevar );
    std::swap( d_data,  rhs.d_data );
    return *this;
}
MathExpr::~MathExpr( )
{
    if ( d_fun != nullptr )
        te_free(d_fun);
}


/************************************************************
*  Initialize the data                                      *
************************************************************/
void MathExpr::initialize( const std::string& expression, const std::vector<std::string>& variables )
{
    d_expr = expression;
    d_vars = variables;
    d_data.resize(variables.size(),0);
    d_tevar.resize(variables.size());
    for (size_t i=0; i<variables.size(); i++)
        d_tevar[i] = { d_vars[i].c_str(), &d_data[i], TE_VARIABLE, nullptr };
    int error = 0;
    d_fun = te_compile( d_expr.c_str(), d_tevar.data(), variables.size(), &error);
    if ( error != 0 ) {
        if ( d_fun != nullptr )
            te_free(d_fun);
        std::string msg = "Error calling te_compile (" + std::to_string(error) + ")" + expression;
        AMP_ERROR(msg);
    }
}


/************************************************************
*  Evaluate the expression                                  *
************************************************************/
double MathExpr::eval( const std::vector<double>& data )
{
    AMP_INSIST(data.size()==d_data.size(),"Number of arguments does not match number of variables");
    for (size_t i=0; i<data.size(); i++)
        d_data[i] = data[i];
    double result = te_eval(d_fun);
    return result;
}


} // AMP namespace

