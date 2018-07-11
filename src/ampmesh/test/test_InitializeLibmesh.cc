#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/initializeLibMesh.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"

#ifdef USE_AMP_VECTORS
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/VectorSelector.h"
#endif

#include "ProfilerApp.h"

#include <sstream>
#include <string>


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
