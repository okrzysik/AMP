#include "AMP/utils/MathExpr.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/tinyexpr/tinyexpr.h"


namespace AMP {


/************************************************************
 *  Constructors/Destructor                                  *
 ************************************************************/
MathExpr::MathExpr() : d_fun( nullptr ) {}
MathExpr::MathExpr( std::string expression, std::vector<std::string> variables )
    : d_expr( std::move( expression ) ), d_fun( nullptr ), d_vars( std::move( variables ) )
{
    initialize();
}
MathExpr::MathExpr( MathExpr &&rhs ) : d_fun( nullptr )
{
    std::swap( d_expr, rhs.d_expr );
    std::swap( d_fun, rhs.d_fun );
    std::swap( d_vars, rhs.d_vars );
    std::swap( d_tevar, rhs.d_tevar );
    std::swap( d_data, rhs.d_data );
}
MathExpr &MathExpr::operator=( MathExpr &&rhs )
{
    if ( this == &rhs )
        return *this;
    std::swap( d_expr, rhs.d_expr );
    std::swap( d_fun, rhs.d_fun );
    std::swap( d_vars, rhs.d_vars );
    std::swap( d_tevar, rhs.d_tevar );
    std::swap( d_data, rhs.d_data );
    return *this;
}
MathExpr::~MathExpr() { te_free( d_fun ); }


/************************************************************
 *  Initialize the data                                      *
 ************************************************************/
void MathExpr::initialize()
{
    d_tevar.resize( d_vars.size() );
    d_data.resize( d_vars.size(), 0 );
    for ( size_t i = 0; i < d_vars.size(); i++ )
        d_tevar[i] = { d_vars[i].c_str(), &d_data[i], TE_VARIABLE, nullptr };
    int error = 0;
    d_fun     = te_compile( d_expr.c_str(), d_tevar.data(), d_vars.size(), &error );
    if ( error != 0 ) {
        te_free( d_fun );
        d_fun           = nullptr;
        std::string msg = "Error calling te_compile (" + std::to_string( error ) + ")" + d_expr;
        AMP_ERROR( msg );
    }
}


/************************************************************
 *  Evaluate the expression                                  *
 ************************************************************/
double MathExpr::operator()( const std::initializer_list<double> &data )
{
    auto it = data.begin();
    for ( size_t i = 0; i < d_data.size(); i++, ++it )
        d_data[i] = *it;
    return te_eval( d_fun );
}
double MathExpr::operator()( const double *data )
{
    for ( size_t i = 0; i < d_data.size(); i++ )
        d_data[i] = data[i];
    return te_eval( d_fun );
}


} // namespace AMP
