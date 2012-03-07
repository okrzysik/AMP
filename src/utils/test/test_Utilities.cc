#include "utils/Utilities.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

#include <string>

#include <assert.h>

#include <fstream>

#include <sys/stat.h>
#include "../Utilities.h"
#include "../AMPManager.h"
#include "../PIO.h"


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
    float mine, wrong, close, relDiff, simpleEps=1e-4;

    // simple.
    mine=1.e-8;
    close=1.000000001e-8;
    wrong=2.e-8;
    
    if(  AMP::Utilities::approx_equal( mine, close, simpleEps ) ) ut->passes("Float precision passes simple check.");
    if( !AMP::Utilities::approx_equal( mine, wrong, simpleEps ) ) ut->passes("Float precision passes simple check.");
    
    mine=1.e-32;
    relDiff = pow( std::numeric_limits<float>::epsilon(), 0.75);
    for (int i=0; i<32; i++) {
      wrong = ( mine ) / ( 1.-relDiff*( 1.01 ) );
      close = ( mine ) / ( 1.-relDiff*( 0.99 ) );
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
    relDiff = pow( std::numeric_limits<float>::epsilon(), 0.75);
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


//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;
    
    try {
        mytest(&ut);
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing test_Utilities, " << err.what() << std::endl;
        ut.failure("test_Utilities");
    } catch( ... ) {
        std::cout << "ERROR: While testing test_Utilities, An unknown exception was thrown." << std::endl;
        ut.failure("test_Utilities");
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   

