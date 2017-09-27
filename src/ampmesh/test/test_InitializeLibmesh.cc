#include <sstream>
#include <string>

#include "ProfilerApp.h"
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/initializeLibMesh.h"
#include "utils/Writer.h"

#ifdef USE_AMP_VECTORS
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#include "vectors/VectorSelector.h"
#endif


int main( int argc, char **argv )
{
    AMP::AMPManager::startup( argc, argv );

    { // Limit scope so variables are destroyed before shutdown
        AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
        AMP::AMP_MPI splitComm = globalComm.split( globalComm.getRank() % 2 );
        AMP_ASSERT( !AMP::Mesh::initializeLibMesh::isInitialized() );
        AMP_ASSERT( AMP::Mesh::initializeLibMesh::canBeInitialized( globalComm ) );
        AMP_ASSERT( AMP::Mesh::initializeLibMesh::canBeInitialized( splitComm ) );
        auto libmesh = AMP::make_shared<AMP::Mesh::initializeLibMesh>( splitComm );
        AMP_ASSERT( AMP::Mesh::initializeLibMesh::isInitialized() );
        if ( globalComm.getSize() > 1 )
            AMP_ASSERT( !AMP::Mesh::initializeLibMesh::canBeInitialized( globalComm ) );
        AMP_ASSERT( AMP::Mesh::initializeLibMesh::canBeInitialized( splitComm ) );
        libmesh.reset();
        AMP_ASSERT( AMP::Mesh::initializeLibMesh::canBeInitialized( globalComm ) );
    }

    AMP::AMPManager::shutdown();
    return 0;
}
