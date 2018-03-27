#ifndef included_AMP_materials_testMaterialHelper
#define included_AMP_materials_testMaterialHelper


#include <string>
#include <vector>

#include "AMP/utils/UnitTest.h"


static const size_t NSUCCESS = 12;
static const size_t NARGEVAL = 2;
static const size_t NVECTOR  = 6;
static const size_t NTENSOR  = 6;


class PropTestResult
{
public:
    PropTestResult()
        : range( false ),
          params( false ),
          unknown( false ),
          name( "none" ),
          isVector( false ),
          isTensor( false )
    {
        for ( auto &elem : success )
            elem = false;
        for ( auto &elem : nargeval )
            elem = false;
        for ( auto &elem : vector )
            elem = false;
        for ( auto &elem : tensor )
            elem = false;
    }
    bool range;
    bool success[NSUCCESS];
    bool params;
    bool nargeval[NARGEVAL];
    bool unknown;
    std::string name;
    bool vector[NVECTOR];
    bool isVector;
    bool tensor[NTENSOR];
    bool isTensor;
};


class MatTestResult
{
public:
    MatTestResult() : undefined( false ), unknown( false ), creationGood( false ), name( "none" ) {}
    void print() const;
    void record( AMP::UnitTest &ut ) const;
    bool undefined;
    bool unknown;
    bool creationGood;
    std::string name;
    std::vector<PropTestResult> propResults;
};


#include "testMaterialHelpers.hpp"

#endif
