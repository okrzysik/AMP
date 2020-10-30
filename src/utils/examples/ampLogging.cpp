#include <iostream>

#include "AMP/utils/AMPManager.h"
#include "AMP/utils/PIO.h"

int main( int argc, char **argv )
{
    // Every AMP program needs to call startup
    AMP::AMPManager::startup( argc, argv );

    // Tells AMP to only log from rank 0 to the file ampLogOnlyNodeZero
    AMP::PIO::logOnlyNodeZero( "ampLogOnlyNodeZero" );

    int rank = AMP::AMP_MPI( AMP_COMM_WORLD ).getRank();

    // after running look in ampLogOnlyNodeZero to see output
    AMP::plog << "Logging from rank " << rank << std::endl;

    // Tells AMP to begin logging from all ranks to the files ampLogAll.*
    AMP::PIO::logAllNodes( "ampLogAll" );

    // after running look in ampAll.* to see output from each rank
    AMP::plog << "Logging from rank " << rank << std::endl;

    // Every AMP program needs to call shutdown
    AMP::AMPManager::shutdown();

    return 0;
}
