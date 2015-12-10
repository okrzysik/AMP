#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <stdexcept>
#include <string.h>
#include <stdint.h>

#include "utils/AMPManager.h"
#include "utils/Utilities.h"
#include "utils/UnitTest.h"
#include "utils/PIO.h"
#include "utils/Array.h"

// The main function
int main( int argc, char *argv[] ) 
{
    // Startup
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    // these are currently not defined for AMP
    //    AMP::Utilities::setAbortBehavior(true,true,true);
    //    AMP::Utilities::setErrorHandlers();

    // Limit the scope of variables
    {
        // Create several matricies
        AMP::Array<double> M1, M2(10,3);
        M1.resize(10,5);
        for (int i=0; i<10; i++) {
            for (int j=0; j<3; j++) {
                M1(i,j) = i+10*j;
                M2(i,j) = i+10*j;
            }
        }
        M1.resize(10,3);
        AMP::Array<double> M3(M1);
        AMP::Array<double> M4 = M2;
        AMP::Array<double> M5 = M1;
        M5(0,0) = -1;
        if ( M1==M2 && M1==M3 && M1==M4 && M1!=M5 ) 
            ut.passes("AMP::Array constructors");
        else
            ut.failure("AMP::Array constructors");
        // Test std::string
        bool pass = true;
        AMP::Array<std::string> S;
        pass = pass && S.length()==0;
        S.resize(1);
        pass = pass && S.length()==1;
        pass = pass && S(0).size()==0;
        S(0) = std::string("test");
        pass = pass && S(0)=="test";
        if ( pass )
            ut.passes("AMP::Array string");
        else
            ut.failure("AMP::Array string");
        // Test a failed allocation
        try {
            size_t N = 10000;
            AMP::Array<double> M(N,N,N);
            #if defined(__APPLE__)
                ut.expected_failure("Failed allocation succeeded (MAC)");
            #else
                ut.failure("Failed allocation succeeded???");
            #endif
            AMP_ASSERT(M.length()==N*N*N);
        } catch (...) {
            ut.passes("Caught failed allocation");
        }
        // Test math opertors
        if ( M1.min()==0 )
            ut.passes("min");
        else
            ut.failure("min");
        if ( M1.max()==29 )
            ut.passes("max");
        else
            ut.failure("max");
        if ( !M1.NaNs() )
            ut.passes("NaNs");
        else
            ut.failure("NaNs");
        // Test find
        std::vector<size_t> index = M1.find( 7, [](double a, double b){return a==b;} );
        if ( index.size()!=1 )
            ut.failure("find");
        else if ( index[0]==7 )
            ut.passes("find");
        else
            ut.failure("find");
        // Test the time required to create a view
        AMP::Array<double> M_view;
        double t1 = AMP::Utilities::time();
        for (size_t i=0; i<100000; i++) {
            M_view.viewRaw( {M1.size(0),M1.size(1)}, M1.data() );
            NULL_USE(M_view);
        }
        double t2 = AMP::Utilities::time();
        if ( M_view == M1 )
            ut.passes("view");
        else
            ut.failure("view");
        AMP::pout << "Time to create view: " << (t2-t1)*1e9/100000 << " ns\n";
    }

    // Finished
    ut.report();
    int num_failed = static_cast<int>(ut.NumFailGlobal());
    if ( num_failed==0 )
        AMP::pout << "All tests passed\n";
    AMP::AMPManager::shutdown();
    return num_failed;
}
