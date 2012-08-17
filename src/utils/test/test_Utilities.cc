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


// This checks approx_equal
template <class T>
void testApproxEqualInt( AMP::UnitTest *ut )
{
    std::string type_name(typeid(T).name());
    if(  AMP::Utilities::approx_equal<T>( 100000, 100000 ) &&
         AMP::Utilities::approx_equal_abs<T>( 100000, 100000 ) && 
         !AMP::Utilities::approx_equal<T>( 100000, 100001 ) &&
         !AMP::Utilities::approx_equal_abs<T>( 100000, 100001 ) ) 
        ut->passes("Integer (" + type_name + ") passes simple check.");
    else
        ut->failure("Integer (" + type_name + ") passes simple check.");
    
    if(  AMP::Utilities::approx_equal_abs<T>( 100001, 100000, 1 ) && 
         !AMP::Utilities::approx_equal_abs<T>( 100002, 100000, 1 ) ) 
        ut->passes("Integer (" + type_name + ") passes close simple check.");
    else
        ut->failure("Integer (" + type_name + ") passes close simple check.");
}
template <class T>
void testApproxEqual( AMP::UnitTest *ut )
{
    std::string type_name(typeid(T).name());

    T mine = 1.0;
    T close_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.8));
    T wrong_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.7));
    T close_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.8);
    T wrong_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.7);
    if( AMP::Utilities::approx_equal( mine, close_rel ) && 
        AMP::Utilities::approx_equal_abs( mine, close_abs ) &&
        !AMP::Utilities::approx_equal( mine, wrong_rel ) && 
        !AMP::Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes(type_name+" passes simple check near 1");
    else
        ut->failure(type_name+" passes simple check near 1");

    mine = 1e-6;
    close_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.8));
    wrong_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.7));
    close_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.8);
    wrong_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.7);
    if( AMP::Utilities::approx_equal( mine, close_rel ) && 
        AMP::Utilities::approx_equal_abs( mine, close_abs ) &&
        !AMP::Utilities::approx_equal( mine, wrong_rel ) && 
        !AMP::Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes(type_name+" passes simple check near 1e-6");
    else
        ut->failure(type_name+" passes simple check near 1e-6");

    mine = -1e-32;
    close_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.8));
    wrong_rel = mine * (1.0+pow(std::numeric_limits<T>::epsilon(),(T)0.7));
    close_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.8);
    wrong_abs = mine + pow(std::numeric_limits<T>::epsilon(),(T)0.7);
    if( AMP::Utilities::approx_equal( mine, close_rel ) && 
        AMP::Utilities::approx_equal_abs( mine, close_abs ) &&
        !AMP::Utilities::approx_equal( mine, wrong_rel ) && 
        !AMP::Utilities::approx_equal_abs( mine, wrong_abs ) )
        ut->passes(type_name+" passes simple check near -1e-32");
    else
        ut->failure(type_name+" passes simple check near -1e-32");

}


// Function to return the call stack
std::vector<std::string> get_call_stack()
{
    return AMP::Utilities::getCallStack();
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
        testApproxEqualInt<int>( &ut );
        testApproxEqualInt<unsigned int>( &ut );
        testApproxEqualInt<size_t>( &ut );
        testApproxEqual<float>( &ut );
        testApproxEqual<double>( &ut );

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
        std::vector<std::string> call_stack = get_call_stack();
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

        // Finished testing, report the results
        ut.report();
        num_failed = ut.NumFailGlobal();
    }

    // Shutdown
    AMP::AMPManager::shutdown();

    // Finished successfully
    return num_failed;
}   

