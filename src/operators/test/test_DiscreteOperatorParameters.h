#ifndef __test_OperatorParameters_h
#define __test_OperatorParameters_h

// Test harness files
#include "utils/UnitTest.h"

#include "test_DiscreteOperatorParameterTests.h"

namespace AMP {
namespace unit_test {

template <typename TEST>
class OperatorParameterTest : public AMP::UnitTest
{
protected:
    explicit OperatorParameterTest( UnitTest *ut ) { d_ut = ut; }
    UnitTest *d_ut;

public:
    virtual ~OperatorParameterTest() {}

    void passes(const std::string &in) override { d_ut->passes(in); }
    void failure(const std::string &in) override { d_ut->failure(in); }
    void expected_failure(const std::string &in) override {
      d_ut->expected_failure(in);
    }

    static void run( UnitTest *ut )
    {
        OperatorParameterTest test( ut );
        TEST::run_test( &test );
    }
};
}
}


#endif
