#ifndef __test_TimeIntegratorParameters_h
#define __test_TimeIntegratorParameters_h

#include "test_TimeIntegratorParameterTests.h"
#include "utils/UnitTest.h"

namespace AMP {
namespace unit_test {

template <typename TEST>
class  TimeIntegratorParameterTest : public AMP::UnitTest
{
protected:
    TimeIntegratorParameterTest ( UnitTest *ut ) { d_ut=ut; }
    UnitTest *d_ut;

public:

    ~TimeIntegratorParameterTest () {}

    void passes(const std::string &in) { d_ut->passes(in); }
    void failure(const std::string &in) { d_ut->failure(in); }
    void expected_failure(const std::string &in) { d_ut->expected_failure(in); }

    static void  run ( UnitTest  *ut )
    {
        try {
            std::cout << "Running " << TEST::get_test_name() << std::endl;
            TimeIntegratorParameterTest test ( ut );
            TEST::run_test ( &test );
        } catch ( std::runtime_error &e ) {
            std::stringstream t;
            t << "Exception thrown: " << e.what();
            ut->failure ( t.str().c_str() );
        } catch ( ... ) {
            ut->failure ( "Unknown exception" );
        }
    }

};

}
}





#endif
