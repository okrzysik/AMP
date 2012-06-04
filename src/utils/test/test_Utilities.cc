#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"


/************************************************************************
*                                                                       *
* This tests whether we can create and use an InputManager object       *
*                                                                       *
************************************************************************/
void mytest(AMP::UnitTest *ut)
{

  // Approx_Equal
  { // double
    double mine, wrong, close, relDiff, simpleEps=1e-4;

    // simple.
    mine=1.e-8;
    close=1.000000001e-8;
    wrong=2.e-8;
    if(  AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->passes("Double precision passes simple check.");
    if( !AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->passes("Double precision passes simple check.");
    
    mine=1.e-32;
    relDiff = pow( std::numeric_limits<double>::epsilon(), 0.75);
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close ) ) { ut->failure("Double precision is checked positively incorrectly for default eps");
      std::cout << i <<"+:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      if(  AMP::Utilities::approx_equal( mine, wrong ) ) { ut->failure("Double precision is checked negatively incorrectly for default eps");
      std::cout << i <<"-:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      mine *=10.;
    }

    mine=1.e-32;
    relDiff=1e-6;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Double precision is checked positively incorrectly for eps = 1e-6");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Double precision is checked negatively incorrectly for eps = 1e-6");
      mine *=10.;
    }

    mine=1.e-32;
    relDiff=1e-14;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) { ut->failure("Double precision is checked positively incorrectly for eps = 1e-14");
      std::cout << i <<"+:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Double precision is checked negatively incorrectly for eps = 1e-14");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff = pow( std::numeric_limits<double>::epsilon(), 0.75);
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close ) ) ut->failure("Negative double precision is checked positively incorrectly for default eps");
      if(  AMP::Utilities::approx_equal( mine, wrong ) ) ut->failure("Negative double precision is checked negatively incorrectly for default eps");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff=1e-6;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Negative double precision is checked positively incorrectly for eps = 1e-6");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Negative double precision is checked negatively incorrectly for eps = 1e-6");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff=1e-14;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Negative double precision is checked positively incorrectly for eps = 1e-14");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Negative double precision is checked negatively incorrectly for eps = 1e-14");
      mine *=10.;
    }
  }
//---------------------------------------------------------------------------//
  { // float
    float mine, wrong, close, relDiff, simpleEps=(float)1e-4;

    // simple.
    mine=(float)1.e-8;
    close=(float)1.000000001e-8;
    wrong=(float)2.e-8;
    
    if(  AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->passes("Float precision passes simple check.");
    if( !AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->passes("Float precision passes simple check.");
    
    mine=(float)1.e-32;
    relDiff = pow( std::numeric_limits<float>::epsilon(), (float) 0.75);
    for (int i=0; i<32; i++) {
      wrong = (float) ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = (float) ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close ) ) { ut->failure("Float precision is checked positively incorrectly for default eps");
      std::cout << i <<"+:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      if(  AMP::Utilities::approx_equal( mine, wrong ) ) { ut->failure("Float precision is checked negatively incorrectly for default eps");
      std::cout << i <<"-:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      mine *=10.;
    }

    mine=1.e-32;
    relDiff=1e-5;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Float precision is checked positively incorrectly for eps = 1e-5");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Float precision is checked negatively incorrectly for eps = 1e-5");
      mine *=10.;
    }

    mine=1.e-32;
    relDiff=1e-3;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) { ut->failure("Float precision is checked positively incorrectly for eps = 1e-3");
      std::cout << i <<"+:"<< mine <<":"<< close <<":"<< wrong <<":"<< relDiff <<std::endl; }
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Float precision is checked negatively incorrectly for eps = 1e-3");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff = pow( std::numeric_limits<float>::epsilon(), (float) 0.75);
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close ) ) ut->failure("Negative float precision is checked positively incorrectly for default eps");
      if(  AMP::Utilities::approx_equal( mine, wrong ) ) ut->failure("Negative float precision is checked negatively incorrectly for default eps");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff=1e-5;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Negative float precision is checked positively incorrectly for eps = 1e-5");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Negative float precision is checked negatively incorrectly for eps = 1e-5");
      mine *=10.;
    }

    mine=-1.e-32;
    relDiff=1e-3;
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
      if( !AMP::Utilities::approx_equal( mine, close, relDiff ) ) ut->failure("Negative float precision is checked positively incorrectly for eps = 1e-3");
      if(  AMP::Utilities::approx_equal( mine, wrong, relDiff ) ) ut->failure("Negative float precision is checked negatively incorrectly for eps = 1e-3");
      mine *=10.;
    }
  }
//---------------------------------------------------------------------------//
  { // integer
    int mine, wrong, close;
    double simpleEps=1e-4;

    // simple.
    mine =100000;
    close=100001;
    wrong=1;
    if(  AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->passes("Integer passes close simple check.");
    if( !AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->passes("Integer passes wrong simple check.");
    
    // zeros
    mine=0;
    close=0;
    wrong=10;
    if( !AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->failure("Integer fails double zero check.");
    if(  AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->failure("Integer fails wrong zero check.");
    
    // their zero
    mine=1;
    wrong=0;
    if(  AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->failure("Integer fails wrong 2nd zero check.");
    
    // same
    mine=1;
    close=1;
    wrong=-1;
    if( !AMP::Utilities::approx_equal( mine, close ) ) ut->failure("Integer fails close 2nd simple check.");
    if(  AMP::Utilities::approx_equal( mine, wrong ) ) ut->failure("Integer fails wrong 2nd simple check.");

    // simple.
    mine =-100000;
    close=-100001;
    wrong= 100000;
    if(  AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->passes("Integer passes close simple check.");
    if( !AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->passes("Integer passes wrong simple check.");
    
    // zeros
    mine=0;
    close=0;
    wrong=-10;
    if( !AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->failure("Integer fails double zero check.");
    if(  AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->failure("Integer fails wrong zero check.");
    
    // their zero
    mine=-1;
    wrong=0;
    if(  AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->failure("Integer fails wrong 2nd zero check.");
    
    // same
    mine=-1;
    close=-1;
    wrong= 1;
    if( !AMP::Utilities::approx_equal( mine, close ) ) ut->failure("Integer fails close 2nd simple check.");
    if(  AMP::Utilities::approx_equal( mine, wrong ) ) ut->failure("Integer fails wrong 2nd simple check.");
 }   
 ut->passes("It reached the end of the test.");
}



//  This test will start and shutdown AMP
int main(int argc, char *argv[])
{

    // Control the behavior of the startup
    AMP::AMPManagerProperties startup_properties;

    // Start AMP
    AMP::AMPManager::startup(argc,argv,startup_properties);
    int num_failed=0;

    // Limit the scope of variables
    { 
        AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

        // Create the unit test
        AMP::UnitTest ut;

        // Try converting an int to a string
        if ( AMP::Utilities::intToString(37,0)=="37" && AMP::Utilities::intToString(37,3)=="037" )
            ut.passes("Convert int to string");
        else
            ut.failure("Convert int to string");

        // Test approx_equal
        if ( AMP::Utilities::approx_equal(1.0,1.0+1e-13) && !AMP::Utilities::approx_equal(1.0,1.0+1e-11) )
            ut.passes("approx_equal");
        else
            ut.failure("approx_equal");

        // Test quicksort performance
        size_t N = 10000;
        std::vector<int> data1(N);
        srand ( time(NULL) );
        for (size_t i=0; i<N; i++)
            data1[i] = rand();
        std::vector<int> data2 = data1;
        std::vector<int> data3 = data1;
        double t1 = AMP::AMP_MPI::time();
        AMP::Utilities::quicksort(data1);
        double t2 = AMP::AMP_MPI::time();
        std::sort(data2.begin(),data2.end());
        double t3 = AMP::AMP_MPI::time();
        std::sort(&data3[0],&data3[0]+data3.size());
        double t4 = AMP::AMP_MPI::time();
        bool pass = true;
        for (size_t i=0; i<N; i++) {
            if ( data1[i]!=data2[i] )
                pass = false;
        }
        if ( pass )
            ut.passes("quicksort sorts correctly");
        else
            ut.failure("quicksort sorts correctly");
        std::cout << "quicksort = " << t2-t1 << ", std::sort = " << t3-t2 << ", std::sort(2) = " << t4-t3 << std::endl;
        
        // Test the hash key
        unsigned int key = AMP::Utilities::hash_char("test");
        if ( key == 2087956275 )
            ut.passes("Got the expected hash key");
        else
            ut.failure("Got the expected hash key");

        // Test the memory usage
        size_t n_bytes = AMP::Utilities::getMemoryUsage();
        if ( globalComm.getRank()==0 )
            std::cout << "Number of bytes used for a basic test: " << n_bytes << std::endl;
        if ( n_bytes > 1e4 )
            ut.passes("getMemoryUsage");
        else
            ut.failure("getMemoryUsage");

        // Test getting the current call stack
        std::vector<std::string> call_stack = AMP::Utilities::getCallStack();
        if ( globalComm.getRank()==0 ) {
            std::cout << "Call stack:" << std::endl;
            for (size_t i=0; i<call_stack.size(); i++)
                std::cout << "   " << call_stack[i];
        }
        if ( call_stack.size() > 0 )
            ut.passes("non empty call stack");
        else
            ut.failure("non empty call stack");

        // Test deleting and checking if a file exists
        if ( globalComm.getRank()==0 ) {
            FILE *fid = fopen( "testDeleteFile.txt", "w" );
            fputs("Temporary test",fid);
            fclose(fid);
            if ( AMP::Utilities::fileExists("testDeleteFile.txt") )
                ut.passes("File exists");
            else
                ut.failure("File exists");
            AMP::Utilities::deleteFile("testDeleteFile.txt");
            if ( !AMP::Utilities::fileExists("testDeleteFile.txt") )
                ut.passes("File deleted");
            else
                ut.failure("File deleted");
        }

		// Run Kevin's test
        try {
            mytest(&ut);
        } catch (std::exception &err) {
            std::cout << "ERROR: While testing test_Utilities, " << err.what() << std::endl;
            ut.failure("test_Utilities");
        } catch( ... ) {
            std::cout << "ERROR: While testing test_Utilities, An unknown exception was thrown." << std::endl;
            ut.failure("test_Utilities");
        }

        // Finished testing, report the results
        ut.report();
        num_failed = ut.NumFailGlobal();
    }

    // Shutdown
    AMP::AMPManager::shutdown();

    // Finished successfully
    return num_failed;
}   

