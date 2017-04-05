#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/MathExpr.h"

#include <math.h>


using namespace AMP;


void test_Expression( const std::string& expression,
    const std::vector<std::string>& vars,
    const std::vector<std::vector<double>>& data,
    const std::vector<double>& ans,
    UnitTest& ut )
{
    bool pass = true;
    MathExpr fun( expression, vars );
    auto fun2 = fun;
    for (size_t i=0; i<ans.size(); i++) {
        auto r1 = fun.eval( data[i] );
        auto r2 = fun.eval( data[i] );
        pass = pass && Utilities::approx_equal_abs( r1, ans[i], 1e-12 );
        pass = pass && Utilities::approx_equal_abs( r1, r2, 1e-14 );
    }
    if ( pass )
        ut.passes(expression);
    else
        ut.failure(expression);
}


int main( int argc, char *argv[] )
{
    AMPManager::startup( argc, argv );
    UnitTest ut;

    // Basic sanity check
    MathExpr fun( "5*8" );
    if ( fun.eval()==40 )
        ut.passes("5*8");
    else
        ut.failure("5*8");

    // Simple check
    fun = MathExpr( "2*x", {"x"} );
    if ( fun.eval( {3.2} ) == 6.4 )
        ut.passes("2*x");
    else
        ut.failure("2*x");

    // Test some basic functions
    test_Expression( "sin(0.8)", {}, {{}}, {sin(0.8)}, ut );
    test_Expression( "3*x", {"x"}, {{3.2}, {6.5}}, {9.6,19.5}, ut );
    test_Expression( "sqrt(x^2+y^2)", {"x","y"}, {{3.0,4.0}, {1.5,3.5}}, {5.0,sqrt(14.5)}, ut );
    test_Expression( "x*exp(y)", {"x","y"}, {{3.0,4.0}, {1.5,3.5}}, {3*exp(4),1.5*exp(3.5)}, ut );
    test_Expression( "a+b+c+d+e+f+g+h", {"a","b","c","d","e","f","g","h"},
        {{1., 2., 3., 4., 5., 6., 7., 8.}}, {36}, ut );

    // Finished
    ut.report();
    int N_errors = ut.NumFailGlobal();
    AMPManager::shutdown();
    return N_errors;
}
