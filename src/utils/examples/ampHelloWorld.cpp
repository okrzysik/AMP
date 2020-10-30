#include "AMP/utils/AMPManager.h"
#include "AMP/utils/PIO.h"

int main( int argc, char **argv )
{
    // Every AMP program starts with a startup call
    // Internally AMP then initializes all the
    // third party libraries it is linking in
    AMP::AMPManager::startup( argc, argv );

    // AMP uses pout to print only from process 0
    // Printing from other ranks is disabled in parallel
    AMP::pout << "Hello World from AMP!!" << std::endl;

    // Every AMP program needs to call shutdown
    // to properly release resources used by AMP
    // and third party libraries it is using
    AMP::AMPManager::shutdown();

    return 0;
}
