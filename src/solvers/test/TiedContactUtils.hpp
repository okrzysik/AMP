
#include <algorithm>
#include <iostream>

void createNodeSet( AMP::Mesh::MeshManager::Adapter::shared_ptr meshAdapter,
                    const unsigned int bndId,
                    std::vector<PointAndId> &nodeSet )
{
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator bnd =
        meshAdapter->beginOwnedBoundary( bndId );
    AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator end_bnd =
        meshAdapter->endOwnedBoundary( bndId );

    nodeSet.clear();

    for ( ; bnd != end_bnd; ++bnd ) {
        PointAndId tmp;
        tmp.d_id    = bnd->globalID();
        tmp.d_pt[0] = bnd->x();
        tmp.d_pt[1] = bnd->y();
        tmp.d_pt[2] = bnd->z();
        nodeSet.push_back( tmp );
    } // end for bnd

    std::sort( nodeSet.begin(), nodeSet.end() );

    std::cout << "Created contact surface nodeset for " << ( meshAdapter->getMeshName() )
              << std::endl;
}

void createMasterSlaveMap( AMP::Mesh::MeshManager::Adapter::shared_ptr masterMeshAdapter,
                           AMP::Mesh::MeshManager::Adapter::shared_ptr slaveMeshAdapter,
                           const unsigned int masterId,
                           const unsigned int slaveId,
                           std::vector<unsigned int> &masterContactNodes,
                           std::vector<unsigned int> &slaveContactNodes,
                           std::vector<unsigned int> &masterVolumeNodes,
                           std::vector<unsigned int> &slaveVolumeNodes )
{
    std::vector<PointAndId> masterContactSet;
    std::vector<PointAndId> slaveContactSet;

    createNodeSet( masterMeshAdapter, masterId, masterContactSet );
    createNodeSet( slaveMeshAdapter, slaveId, slaveContactSet );

    std::cout << "Master contact surface has " << ( masterContactSet.size() ) << " nodes."
              << std::endl;
    std::cout << "Slave contact surface has " << ( slaveContactSet.size() ) << " nodes."
              << std::endl;

    const size_t numContactNodes = masterContactSet.size();

    AMP_INSIST( ( masterContactSet.size() == slaveContactSet.size() ),
                "Masterset size differs from Slaveset size" );

    for ( size_t i = 0; i < numContactNodes; i++ ) {
        double mX        = masterContactSet[i].d_pt[0];
        double mY        = masterContactSet[i].d_pt[1];
        double mZ        = masterContactSet[i].d_pt[2];
        double sX        = slaveContactSet[i].d_pt[0];
        double sY        = slaveContactSet[i].d_pt[1];
        double sZ        = slaveContactSet[i].d_pt[2];
        unsigned int mId = masterContactSet[i].d_id;
        unsigned int sId = slaveContactSet[i].d_id;
        char msg[300];
        sprintf( msg,
                 "Master point (%g, %g, %g) and Slave point (%g, %g, %g) differ.\n Master ID = %u "
                 "and Slave ID = %u.\n",
                 mX,
                 mY,
                 mZ,
                 sX,
                 sY,
                 sZ,
                 mId,
                 sId );
        AMP_INSIST( ( masterContactSet[i].pointEquals( slaveContactSet[i] ) ), msg );
    }

    masterContactNodes.resize( numContactNodes );
    slaveContactNodes.resize( numContactNodes );

    for ( size_t i = 0; i < numContactNodes; i++ ) {
        masterContactNodes[i] = masterContactSet[i].d_id;
    }

    for ( size_t i = 0; i < numContactNodes; i++ ) {
        slaveContactNodes[i] = slaveContactSet[i].d_id;
    }

    std::cout << "Master Volume Node Size = ";
    std::cout << ( masterMeshAdapter->numLocalNodes() ) << std::endl;

    std::cout << "Slave Volume Node Size = ";
    std::cout << ( slaveMeshAdapter->numLocalNodes() ) << std::endl;

    masterVolumeNodes.resize( masterMeshAdapter->numLocalNodes() );
    slaveVolumeNodes.resize( slaveMeshAdapter->numLocalNodes() );

    for ( size_t i = 0; i < masterVolumeNodes.size(); i++ ) {
        masterVolumeNodes[i] = i;
    }

    for ( size_t i = 0; i < slaveVolumeNodes.size(); i++ ) {
        slaveVolumeNodes[i] = i;
    }

    for ( size_t i = 0; i < numContactNodes; i++ ) {
        AMP_ASSERT( masterContactNodes[i] < masterVolumeNodes.size() );
        AMP_ASSERT( slaveContactNodes[i] < slaveVolumeNodes.size() );
        masterVolumeNodes[masterContactNodes[i]] = static_cast<unsigned int>( -1 );
        slaveVolumeNodes[slaveContactNodes[i]]   = static_cast<unsigned int>( -1 );
    }

    std::vector<unsigned int> tmpMaster;
    for ( size_t i = 0; i < masterVolumeNodes.size(); i++ ) {
        if ( masterVolumeNodes[i] != ( static_cast<unsigned int>( -1 ) ) ) {
            tmpMaster.push_back( masterVolumeNodes[i] );
        }
    }

    std::vector<unsigned int> tmpSlave;
    for ( size_t i = 0; i < slaveVolumeNodes.size(); i++ ) {
        if ( slaveVolumeNodes[i] != ( static_cast<unsigned int>( -1 ) ) ) {
            tmpSlave.push_back( slaveVolumeNodes[i] );
        }
    }

    masterVolumeNodes = tmpMaster;
    slaveVolumeNodes  = tmpSlave;

    std::cout << "Created Master-Slave map for Tied-Contact." << std::endl;
}
