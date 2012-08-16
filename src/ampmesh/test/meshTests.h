#ifndef included_MeshTests
#define included_MeshTests

#include <vector>
#include <set>

#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MultiMesh.h"
#include "ampmesh/SubsetMesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/MeshIterator.h"
#include "ampmesh/MeshElementVectorIterator.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
#endif


// Helper function to create a map from the base mesh communicator rank to the main mesh communicator
std::map<AMP::Mesh::MeshID,std::vector<int> >  createRankMap( AMP::Mesh::Mesh::shared_ptr mesh )
{
    std::map<AMP::Mesh::MeshID,std::vector<int> > proc_map;
    std::vector<AMP::Mesh::MeshID> meshIDs = mesh->getBaseMeshIDs();
    std::vector<std::pair<int,int> > tmp(mesh->getComm().getSize());
    for (size_t i=0; i<meshIDs.size(); i++) {
        AMP::Mesh::Mesh::shared_ptr mesh2 = mesh->Subset(meshIDs[i]);
        int N_send = 0;
        std::pair<int,int> map;
        if ( mesh2.get() != NULL ) {
            map = std::pair<int,int>(mesh2->getComm().getRank(),mesh->getComm().getRank());
            N_send = 1;
        }
        int N = mesh->getComm().allGather( &map, N_send, &tmp[0] );
        std::vector<int> rank(N);
        for (int j=0; j<N; j++)
            rank[tmp[j].first] = tmp[j].second;
        proc_map.insert( std::pair<AMP::Mesh::MeshID,std::vector<int> >( meshIDs[i], rank ) );
    }
    return proc_map;
}


// This test checks a single mesh element iterator
// ut           Unit test class to report the results
// iterator     local iterator over elements
// N_local      number of local elements for the iterator
// N_ghost      number of ghost elements for the iterator
void ElementIteratorTest( AMP::UnitTest *ut, AMP::Mesh::Mesh::shared_ptr mesh, AMP::Mesh::MeshIterator iterator, 
    const size_t N_local, const size_t N_ghost, const AMP::Mesh::GeomType type )
{
    // For each mesh, get a mapping of it's processor id's to the comm of the mesh
    std::map<AMP::Mesh::MeshID,std::vector<int> > proc_map = createRankMap( mesh );
    // Check that we can get the begin and end iterator
    AMP::Mesh::MeshIterator  begin_it = iterator.begin();
    AMP::Mesh::MeshIterator  end_it = iterator.end();
    if ( N_local+N_ghost==0 ) {
        if ( begin_it == end_it )
            ut->passes("trival iterator begin and end returned");
        else
            ut->failure("trival iterator begin and end returned");
        return;
    }
    if ( begin_it != end_it )
        ut->passes("iterator begin and end returned");
    else
        ut->failure("iterator begin and end returned");

    // Check that the iterator iterates through the proper number of elements
    if ( iterator.size() == N_local+N_ghost )
        ut->passes("regular iterator size()");
    else
        ut->failure("regular iterator size()");
    size_t number_of_local_elements = 0;
    size_t number_of_ghost_elements = 0;
    std::set<AMP::Mesh::MeshElementID>  ids;
    AMP::Mesh::MeshIterator  cur_it = iterator.begin();
    while ( cur_it != end_it ) {
        AMP::Mesh::MeshElementID id = cur_it->globalID();
        ids.insert ( id );
        if ( id.is_local() )
            number_of_local_elements++;
        else
            number_of_ghost_elements++;
        ++cur_it;   // Pre-increment is faster than post-increment
    }
    if ( number_of_local_elements==N_local && number_of_ghost_elements==N_ghost )
        ut->passes("regular iterator count");
    else
        ut->failure("regular iterator count");
    if ( ids.size() == N_local+N_ghost )
        ut->passes("regular iterator uniqueness");
    else
        ut->failure("regular iterator uniqueness");

    // Check that we can increment and decrement properly
    if ( iterator.size()>=2 ) {
        bool pass = true;
        AMP::Mesh::MeshIterator it1 = iterator.begin();
        AMP::Mesh::MeshIterator it2 = iterator.begin();
        AMP::Mesh::MeshIterator it3 = iterator.begin();
        it1++;
        ++it1;
        it2 = it2+2;
        it3+=2;
        if ( it1!=it2 || it1!=it3 )
            pass = false;
        /*it1--;
        --it1;
        it2 = it2-2;
        it3-=2;
        if ( it1!=iterator.begin() || it2!=iterator.begin() || it3!=iterator.begin() )
            pass = false;*/
        if ( pass )
            ut->passes("regular iterator increments/decrements");
        else
            ut->failure("regular iterator increments/decrements");
    }

    // Run element tests
    bool id_pass = true;
    bool type_pass = true;
    bool volume_pass = true;
    bool coord_pass = true;
    bool centroid_pass = true;
    bool elements_pass = true;
    int neighbor_pass = 1;
    cur_it = iterator.begin();
    int myRank = mesh->getComm().getRank();
    int maxRank = mesh->getComm().getSize()-1;
    while ( cur_it != end_it ) {
        AMP::Mesh::MeshElement element = *cur_it;
        // Get the current id
        AMP::Mesh::MeshElementID id = element.globalID();
        if ( id != cur_it->globalID() )
            id_pass = false;
        // Get the owner rank
        AMP::Mesh::MeshID meshID = id.meshID();
        const std::vector<int> &map = proc_map.find(meshID)->second;
        int ownerRank = map[id.owner_rank()];
        // Perform some simple checks
        if ( element.elementType() != type )
            type_pass = false;
        if ( type==AMP::Mesh::Vertex ) {
            std::vector<double> coord = element.coord();
            if ( coord.size()!=mesh->getDim() )
                coord_pass = false;
        } else {
            if ( element.volume() <= 0.0 )
                volume_pass = false;
        }
        std::vector<double> centroid = element.centroid();
        if ( centroid.size()!=mesh->getDim() )
            centroid_pass = false;
        if ( type==AMP::Mesh::Vertex ) {
            std::vector<double> coord = element.coord();
            for (size_t i=0; i<centroid.size(); i++) {
                if ( centroid[i]!=coord[i] ) 
                    centroid_pass = false;
            }
        }
        if ( id.is_local() ) {
            for (int i=0; i<=(int)type; i++) {
                AMP::Mesh::GeomType type2 = (AMP::Mesh::GeomType) i;
                std::vector<AMP::Mesh::MeshElement> pieces = element.getElements(type2);
                if ( pieces.empty() )
                    elements_pass = false;  // Did not return anything
                AMP::Utilities::quicksort(pieces);
                for (size_t j=1; j<pieces.size(); j++) {
                    if ( pieces[j]==pieces[j-1] )
                        elements_pass = false;  // Repeated elements
                }
            }
            std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = element.getNeighbors();
            if ( neighbors.empty() ) {
                if ( element.elementType()==AMP::Mesh::Vertex || element.elementType()==mesh->getDim() )
                    neighbor_pass = 0;
                else if ( neighbor_pass==1 )
                    neighbor_pass = 2;
            }
            if ( ownerRank!=myRank )
                id_pass = false;
        } else {
            if ( ownerRank>maxRank || ownerRank<0 || ownerRank==myRank )
                id_pass = false;
        }
        ++cur_it;   // Pre-increment is faster than post-increment
    }
    if ( id_pass && type_pass && volume_pass && coord_pass && elements_pass && neighbor_pass==1 ) {
        ut->passes( "elements passed" );
    } else {
        if ( !id_pass )
            ut->failure( "elements failed id test" );
        if ( !type_pass )
            ut->failure( "elements failed type test" );
        if ( !volume_pass )
            ut->failure( "elements failed volume test" );
        if ( !coord_pass )
            ut->failure( "elements failed coord test" );
        if ( !centroid_pass )
            ut->failure( "elements failed centroid test" );
        if ( !elements_pass )
            ut->failure( "elements failed getElements test" );
        if ( neighbor_pass==0 ) 
            ut->failure( "elements failed getNeighbors test" );
        else if ( neighbor_pass==2 ) 
            ut->expected_failure( "elements failed getNeighbors test" );
    }
    // Check that we can get the element from the global id for all elements
    cur_it = iterator.begin();
    bool getElem_pass = true;
    for (size_t i=0; i<cur_it.size(); i++) {
        AMP::Mesh::MeshElementID id1 = cur_it->globalID();
        AMP::Mesh::MeshElement elem = mesh->getElement(id1);
        AMP::Mesh::MeshElementID id2 = elem.globalID();
        if ( id1 != id2 )
            getElem_pass = false;
        ++cur_it;
    }
    if ( getElem_pass )
        ut->passes( "Got elements from element ids" );
    else
        ut->failure( "Got elements from element ids" );
}


// Check the different mesh element iterators
void MeshIteratorTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    char message[1000];
    // Loop through different ghost widths
    for (int gcw=0; gcw<=1; gcw++) {
        // Loop through the different geometric entities
        for (int i=0; i<=(int)mesh->getGeomType(); i++) {
            AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) i;
            // Try to create the iterator
            size_t N_local = 0;
            size_t N_ghost = 0;
            AMP::Mesh::MeshIterator iterator;
            bool iterator_created = true;
            try {
                N_local = mesh->numLocalElements(type);
                N_ghost = mesh->numGhostElements(type,gcw);
                iterator = mesh->getIterator(type,gcw);
                sprintf(message,"Element iterator created (gcw=%i)",gcw);
                ut->passes(message);
            } catch (...) {
                iterator_created = false;
                if ( i==0 ) {
                    sprintf(message,"Node iterator failed (gcw=%i)",gcw);
                    ut->failure(message);
                } else if ( type==mesh->getGeomType() ) {
                    sprintf(message,"Geometric element iterator failed (gcw=%i)",gcw);
                    ut->failure(message);
                } else {
                    sprintf(message,"Intermediate element iterator failed (gcw=%i)",gcw);
                    ut->expected_failure(message);
                }
            }
            // Test the regular iterator over local elements
            if ( iterator_created )
                ElementIteratorTest( ut, mesh, iterator, N_local, N_ghost, type );
            // Add tests with gcw != 0
            // Add const iterator tests
        }
    }    
}


// Test operator operations for iterator
void MeshIteratorOperationTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    // Create some iterators to work with
    AMP::Mesh::MeshIterator A = mesh->getIterator(AMP::Mesh::Vertex,1);
    AMP::Mesh::MeshIterator B = mesh->getIterator(mesh->getGeomType(),0);
    boost::shared_ptr<std::vector<AMP::Mesh::MeshElement> > elements( 
        new std::vector<AMP::Mesh::MeshElement>(A.size()) );
    AMP::Mesh::MeshIterator tmp = A.begin();
    for (size_t i=0; i<A.size(); i++) {
        (*elements)[i] = *tmp;
        ++tmp;
    }
    AMP::Mesh::MeshIterator C = AMP::Mesh::MultiVectorIterator( elements );

    // Check operator== and operator!=
    if ( A==A && B==B && C==C )
        ut->passes("Iterator == with same iterator");
    else
        ut->failure("Iterator == with same iterator");
    if ( !(A!=A) && !(B!=B) && !(C!=C) )
        ut->passes("Iterator != with same iterator");
    else
        ut->failure("Iterator != with same iterator");
    if ( !(A==B) )
        ut->passes("Iterator == with same type, different iterator");
    else
        ut->failure("Iterator == with same type, different iterator");
    if ( A!=B )
        ut->passes("Iterator != with same type, different iterator");
    else
        ut->failure("Iterator != with same type, different iterator");
    if ( A==C && C==A && !(B==C) && !(C==B) )
        ut->passes("Iterator == with different type");
    else
        ut->failure("Iterator == with different type");
    if ( !(A!=C) && !(C!=A) && B!=C && C!=B )
        ut->passes("Iterator != with different type");
    else
        ut->failure("Iterator != with different type");
}


// Test set operations for the iterators
void MeshIteratorSetOPTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    AMP::Mesh::MeshIterator A = mesh->getIterator(AMP::Mesh::Vertex,1);
    AMP::Mesh::MeshIterator B = mesh->getIterator(AMP::Mesh::Vertex,0);
    AMP::Mesh::MeshIterator C = AMP::Mesh::MeshIterator();
    AMP::Mesh::MeshIterator R1, R2;
    // Check Union
    R1 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, A, B);
    R2 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Union, B, C);
    if ( R1.size()==A.size() && R2.size()==B.size() )
        ut->passes("Union iterator create");
    else
        ut->failure("Union iterator create");
    // Check Intersection
    R1 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Intersection, A, B);
    R2 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Intersection, B, C);
    if ( R1.size()==B.size() && R2.size()==0 )
        ut->passes("Intersection iterator create");
    else
        ut->failure("Intersection iterator create");
    // Check Complement
    R1 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Complement, A, B);
    R2 = AMP::Mesh::Mesh::getIterator( AMP::Mesh::Complement, B, C);
    if ( R1.size()==mesh->numGhostElements(AMP::Mesh::Vertex,1) && R2.size()==B.size() )
        ut->passes("Complement iterator create");
    else
        ut->failure("Complement iterator create");
}


// Test the number of elements in the mesh
void MeshCountTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    AMP::AMP_MPI comm = mesh->getComm();
    for (int i=0; i<=(int)mesh->getGeomType(); i++) {
        AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) i;
        size_t N_local = mesh->numLocalElements(type);
        size_t N_global = mesh->numGlobalElements(type);
        size_t N_ghost0 = mesh->numGhostElements(type,0);
        size_t N_ghost1 = comm.sumReduce( mesh->numGhostElements(type,1) );
        size_t N_sum = comm.sumReduce(N_local);
        if ( N_global > 0 )
            ut->passes("Non-trival mesh created");
        else
            ut->failure("Non-trival mesh created");
        if ( N_sum == N_global )
            ut->passes("Sum of local mesh counts matches global count");
        else
            ut->failure("Sum of local mesh counts matches global count");
        if ( N_ghost0 == 0 )
            ut->passes("gcw=0 has no ghost elements");
        else
            ut->failure("gcw=0 has no ghost elements");
        std::vector<AMP::Mesh::MeshID> ids = mesh->getBaseMeshIDs();
        bool is_base_mesh = ids.size()==1 && ids[0]==mesh->meshID();
        if ( N_local != N_global && is_base_mesh ) {
            if ( N_ghost1 > 0 )
                ut->passes("gcw=1 has ghost elements");
            else
                ut->failure("gcw=1 has ghost elements");
        }
    }
}


// Test some basic Mesh properties
void MeshBasicTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    // test that we can get the mesh ID
    AMP::Mesh::MeshID meshID = mesh->meshID();
    if ( meshID>0 && meshID!=AMP::Mesh::MeshID() ) 
        ut->passes("got meshID");
    else
        ut->failure("got meshID");
    // Test that we can subset the mesh for it's self using the meshID
    AMP::Mesh::Mesh::shared_ptr mesh2 = mesh->Subset(meshID);
    if ( mesh2.get()==mesh.get() )
        ut->passes("subset on meshID for self");
    else
        ut->failure("subset on meshID for self");
    // test that we can get and set the mesh name
    std::string meshName = mesh->getName();
    mesh->setName("testing mesh name");
    bool setName = mesh->getName().compare("testing mesh name")==0;
    mesh->setName(meshName);
    if ( meshName.compare("NULL")!=0 ) 
        ut->passes("non-null mesh name");
    else
        ut->failure("non-null mesh name");
    if ( setName ) 
        ut->passes("get/set mesh name");
    else
        ut->failure("get/set mesh name");
    // Test that we can subset the mesh by the mesh name
    mesh2 = mesh->Subset(meshName);
    if ( mesh2.get()==mesh.get() )
        ut->passes("subset on mesh name for self");
    else
        ut->failure("subset on mesh name for self");
    mesh2 = mesh->Subset("Garbage name");
    if ( mesh2.get()==NULL )
        ut->passes("subset on mesh name for garbage");
    else
        ut->failure("subset on mesh name for garbage");
    // Check that the bounding box matches on all processors
    std::vector<double> box1 = mesh->getBoundingBox();
    std::vector<double> box2 = box1;
    mesh->getComm().bcast( &box2[0], (int) box1.size(), 0 );
    bool box_match = true;
    for (size_t i=0; i<box1.size(); i++) {
        if ( box1[i]!=box2[i] ) 
            box_match = false;
    }
    if ( box_match )
        ut->passes("mesh->getBoundingBox returns global bounding box");
    else
        ut->failure("mesh->getBoundingBox returns global bounding box");
}


// This tests checks that all ghost elements are owned by "owner processor"
void VerifyGhostIsOwned( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    for (int type=0; type<=(int)mesh->getGeomType(); type++) {
        int gcw = mesh->getMaxGhostWidth();
        // Build a list of the owned and ghost elements
        std::vector<AMP::Mesh::MeshElementID> owned, ghost;
        owned.reserve( mesh->numLocalElements( (AMP::Mesh::GeomType) type ) );
        ghost.reserve( mesh->numGhostElements( (AMP::Mesh::GeomType) type, gcw ) );
        AMP::Mesh::MeshIterator iterator = mesh->getIterator( (AMP::Mesh::GeomType) type, gcw );
        for (size_t i=0; i<iterator.size(); i++) {
            AMP::Mesh::MeshElementID id = iterator->globalID();
            if ( id.is_local() )
                owned.push_back(id);
            else
                ghost.push_back(id);
            ++iterator;
        }
        // Broadcast the list of ghost ids to everybody
        size_t N_ghost_global = mesh->getComm().sumReduce( ghost.size() );
        if ( N_ghost_global==0 )
            continue;
        std::vector<AMP::Mesh::MeshElementID> ghost_global(N_ghost_global);
        AMP::Mesh::MeshElementID *send_data=NULL;
        if ( ghost.size() > 0 ) { send_data = &ghost[0]; }
        AMP::Mesh::MeshElementID *recv_data = &ghost_global[0];
        mesh->getComm().allGather( send_data, (int) ghost.size(), recv_data );
        // Check that each ghost appears in the owner's rank's list
        AMP::Utilities::quicksort( owned );         // Sort for search
        AMP::Utilities::quicksort( ghost_global );  // Sort for speed
        std::vector<int> found(ghost_global.size(),0);
        AMP::Mesh::Mesh::shared_ptr my_mesh = mesh;
        unsigned int my_rank = my_mesh->getComm().getRank();
        AMP::Mesh::MeshID my_mesh_id = my_mesh->meshID();
        for (size_t i=0; i<N_ghost_global; i++) {
            // Get the current mesh
            if ( ghost_global[i].meshID() != my_mesh_id ) {
                my_mesh_id = ghost_global[i].meshID();
                my_mesh = mesh->Subset( my_mesh_id );
                if ( my_mesh.get()==NULL) 
                    continue;
                my_rank = my_mesh->getComm().getRank();
            }
            // Check if we are the owning rank
            if ( my_mesh.get()==NULL) 
                continue;
            if ( ghost_global[i].owner_rank() != my_rank )
                continue;
            // Check if we have the element
            size_t index = AMP::Utilities::findfirst( owned, ghost_global[i] );
            if ( index == owned.size() ) { index--; }
            if ( owned[index] == ghost_global[i] )
                found[i] = 1;
        }
        mesh->getComm().maxReduce( &found[0], (int) found.size() );
        bool all_found = true;
        for (size_t i=0; i<found.size(); i++) {
            if ( found[i]==0 )
                all_found = false;
        }
        if ( all_found )
            utils->passes("All ghosts are owned by somebody");
        else
            utils->failure("All ghosts are owned by somebody");
    }
}


// This tests loops over all boundary ids
void VerifyBoundaryIDNodeIterator( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    const std::vector<int> bids = mesh->getBoundaryIDs();
    for (size_t i=0; i<bids.size(); i++) {
        int bid = bids[i];
        for (int gcw=0; gcw<=0; gcw++) {
            // Get the iterator over the current boundary id
            AMP::Mesh::MeshIterator curNode = mesh->getBoundaryIDIterator( AMP::Mesh::Vertex, bid, gcw );
            AMP::Mesh::MeshIterator endNode = curNode.end();
            // Get the set of all nodes in the iterator
            bool testPassed = true;
            std::set<AMP::Mesh::MeshElementID>  node_ids;
            while ( curNode != endNode ) {
                node_ids.insert( curNode->globalID() );
                if ( !curNode->isOnBoundary(bid) )
                    testPassed = false;
                curNode++;
            }
            size_t total_size = mesh->getComm().sumReduce(node_ids.size());
            if (total_size==0 )
                testPassed = false;
            // Verify that all nodes were found
            size_t  numFound = 0;
            AMP::Mesh::MeshIterator  curMNode = mesh->getIterator(AMP::Mesh::Vertex,gcw);
            AMP::Mesh::MeshIterator  endMNode = curMNode.end();
            while ( curMNode != endMNode ) {
                if ( curMNode->isOnBoundary( bid ) ) {
                    numFound++;
                    if ( node_ids.find( curMNode->globalID() ) == node_ids.end() )
                        testPassed = false;
                }
                curMNode++;
            }
            if ( numFound != node_ids.size() )
                testPassed = false;
            if ( testPassed )
                utils->passes ( "Found all boundary nodes" );
            else
                utils->failure( "Found all boundary nodes" );
        }
    }
}


// This tests loops over the boundary
void VerifyBoundaryIterator( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    for (int gcw=0; gcw<=0; gcw++) {
        for (int type2=0; type2<=(int)mesh->getGeomType(); type2++) {
            AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) type2;
            // Get the iterator over the current boundary id
            AMP::Mesh::MeshIterator iterator = mesh->getSurfaceIterator( type, gcw );
            size_t global_size = mesh->getComm().sumReduce(iterator.size());
            bool passes = global_size>0;
            if ( boost::dynamic_pointer_cast<AMP::Mesh::SubsetMesh>(mesh).get()==NULL ) {
                if ( mesh->numGlobalElements(type) >= 100 )
                    passes = passes && global_size<mesh->numGlobalElements(type);
            }
            if ( passes )
                utils->passes("Non-trivial surface iterator created");
            else
                utils->failure("Non-trivial surface iterator created");
        }
    }
}


// This tests checks the block ids
void testBlockIDs( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    const std::vector<int> blockIDs = mesh->getBlockIDs();
    if ( blockIDs.size() > 0 )
        utils->passes("Block ids found");
    else if ( (int)mesh->getGeomType() != mesh->getDim() ) 
        utils->expected_failure("Block ids need work for surface meshes");
    else
        utils->failure("Block ids found");
}


// This tests basic id info
void testID( AMP::UnitTest *utils )
{
    unsigned int num_failed0 = utils->NumFailLocal();
    // Create some IDs for testing
    AMP::Mesh::MeshElementID id0;
    AMP::Mesh::MeshElementID id1(false,AMP::Mesh::Vertex,2,1,103);
    AMP::Mesh::MeshElementID id2(true,AMP::Mesh::Vertex,2,1,103);
    AMP::Mesh::MeshElementID id3(true,AMP::Mesh::Volume,2,1,103);
    AMP::Mesh::MeshElementID id4(true,AMP::Mesh::Vertex,3,1,103);
    AMP::Mesh::MeshElementID id5(true,AMP::Mesh::Vertex,2,4,103);
    AMP::Mesh::MeshElementID id6(true,AMP::Mesh::Vertex,2,1,105);
    // Test the default values
    if ( id0.meshID()!=0xFFFFFFFFFFFFFFFF || id0.is_local() || id0.type()!=AMP::Mesh::null || 
        id0.owner_rank()!=0 || id0.local_id()!=0xFFFFFFFF )
        utils->failure("MeshElementID test defaults");
    // Test == and != operators
    if ( !(id1==id1) || !(id1==id2) )
        utils->failure("MeshElementID test ==");
    if ( (id1!=id1) || (id1!=id2) )
        utils->failure("MeshElementID test !=");
    if ( id1==id3 || id1==id4 || id1==id5 || id1==id6 )
        utils->failure("MeshElementID test == (2)");
    // Test that the basic properties were assigned correctly
    if ( id1.is_local() || !id2.is_local() )
        utils->failure("MeshElementID test is_local");
    if ( id1.type()!=AMP::Mesh::Vertex || id1.local_id()!=2 || id1.owner_rank()!=1 || id1.meshID()!=103 )
        utils->failure("MeshElementID test values");
    id1.set_is_local(true);
    id2.set_is_local(false);
    if ( !id1.is_local() || id2.is_local() )
        utils->failure("MeshElementID test is_local (2)");
    // test greater than and less than operators
    if ( !(id1<=id2) || !(id1>=id2) || id1>id2 || id1<id2 )
        utils->failure("MeshElementID test <,> (1)");
    if ( id1>id3 || id1>=id3 || id3<id1 || id3<=id1 )
        utils->failure("MeshElementID test <,> (2)");
    if ( id1>id4 || id1>=id4 || id4<id1 || id4<=id1 )
        utils->failure("MeshElementID test <,> (3)");
    if ( id1>id5 || id1>=id5 || id5<id1 || id5<=id1 )
        utils->failure("MeshElementID test <,> (4)");
    if ( id1>id6 || id1>=id6 || id6<id1 || id6<=id1 )
        utils->failure("MeshElementID test <,> (5)");
    // The elements should sort by meshID, processor id, type, then local id
    std::vector<AMP::Mesh::MeshElementID> list(6);
    list[0] = id0;
    list[1] = id3;
    list[2] = id4;
    list[3] = id6;
    list[4] = id5;
    list[5] = id1;
    AMP::Utilities::quicksort(list);
    if ( list[0]!=id1 || list[1]!=id4 || list[2]!=id3 || list[3]!=id5 || list[4]!=id6 || list[5]!=id0 )
        utils->failure("MeshElementID test sort");
    if ( num_failed0 == utils->NumFailLocal() )
        utils->passes("MeshElementID tests");
    else
        utils->failure("MeshElementID tests");
}


// Test if we correctly identify the node neighbors
void getNodeNeighbors( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    std::map< AMP::Mesh::MeshElementID, std::vector<AMP::Mesh::MeshElementID> > neighbor_list;
    // Get a list of all neighors for each local node
    AMP::Mesh::MeshIterator nodeIterator = mesh->getIterator(AMP::Mesh::Vertex,0);
    std::vector<AMP::Mesh::MeshElementID> neighbors(100);
    for (size_t i=0; i<nodeIterator.size(); i++) {
        std::vector<AMP::Mesh::MeshElement::shared_ptr> elements = nodeIterator->getNeighbors();
        // Store the neighbor list
        neighbors.resize(0);
        for (size_t j=0; j<elements.size(); j++) {
            if ( elements[j].get() != NULL )
                neighbors.push_back(elements[j]->globalID());
        }
        // Sort the neighbor list for easy searching
        AMP::Utilities::quicksort(neighbors);
        std::pair< AMP::Mesh::MeshElementID, std::vector<AMP::Mesh::MeshElementID> > entry;
        entry.first = nodeIterator->globalID();
        entry.second = neighbors;
        neighbor_list.insert( entry );
        ++nodeIterator;
    }
    // First check if the neighbor lists are unique and don't contain self
    {
        std::map< AMP::Mesh::MeshElementID, std::vector<AMP::Mesh::MeshElementID> >::iterator iterator;
        bool contains_self = false;
        bool contains_duplicate = false;
        for (iterator=neighbor_list.begin(); iterator!=neighbor_list.end(); iterator++) {
            std::vector<AMP::Mesh::MeshElementID> neighbors = iterator->second;
            for (size_t i=0; i<neighbors.size(); i++) {
                if ( neighbors[i] == iterator->first )
                    contains_self = true;
                for (size_t j=0; j<i; j++) {
                    if ( neighbors[j] == neighbors[i] )
                        contains_duplicate = true;
                }
            }
        }
        if ( !contains_self )
            utils->passes("Neighbor nodes does not contain self");
        else
            utils->failure("Neighbor nodes does not contain self");
        if ( !contains_duplicate )
            utils->passes("Neighbor nodes does not contain duplicates");
        else
            utils->failure("Neighbor nodes does not contain duplicates");
    }
    // If there are ghost nodes, then some of them must be neighbors
    if ( mesh->numGhostElements(AMP::Mesh::Vertex,1) > 0 ) {
        std::map< AMP::Mesh::MeshElementID, std::vector<AMP::Mesh::MeshElementID> >::iterator iterator;
        bool ghost_neighbors = false;
        for (iterator=neighbor_list.begin(); iterator!=neighbor_list.end(); iterator++) {
            std::vector<AMP::Mesh::MeshElementID> neighbors = iterator->second;
            for (size_t i=0; i<neighbors.size(); i++) {
                if ( !neighbors[i].is_local() )
                    ghost_neighbors = true;
            }
        }
        if ( ghost_neighbors )
            utils->passes("Found ghost neighbor nodes");
        else
            utils->failure("Found ghost neighbor nodes");
    }
    // Loop through all elements (including ghosts), for each owned node, check that all other nodes are neighbors
    bool passed = true;
    AMP::Mesh::MeshIterator elementIterator = mesh->getIterator(mesh->getGeomType(),1);
    for (size_t i=0; i<elementIterator.size(); i++) {
        std::vector<AMP::Mesh::MeshElement> nodes = elementIterator->getElements(AMP::Mesh::Vertex);
        for (size_t j=0; j<nodes.size(); j++) {
            if ( !nodes[j].globalID().is_local() )
                continue;   // Node is not owned, move on
            std::map< AMP::Mesh::MeshElementID, std::vector<AMP::Mesh::MeshElementID> >::iterator iterator;
            iterator = neighbor_list.find( nodes[j].globalID() );
            if ( iterator==neighbor_list.end() ) {
                passed = false;
                break;
            }
            const std::vector<AMP::Mesh::MeshElementID> &neighbors = iterator->second;
            if ( neighbors.size()==0 ) {
                passed = false;
                break;
            }
            for (size_t k=0; k<nodes.size(); k++) {
                if ( k==j )
                    continue;
                size_t index = AMP::Utilities::findfirst( neighbors, nodes[k].globalID() );
                if ( index==neighbors.size() )
                    passed = false;
            }
        }
        ++elementIterator;
    }
    if ( passed )
        utils->passes("Node neighbors found all neighbors");
    else
        utils->failure("Node neighbors found all neighbors");
    // Transfer all neighbor lists to all processors and check that every neighbor node 
    // also has the current node as a neighbor (not finished)
}


// Test the displacement of the mesh
void DisplaceMesh( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    // Test the scalar displacement
    std::vector<double> box1 = mesh->getBoundingBox();
    std::vector<double> displacement(mesh->getDim(),1.0);
    mesh->displaceMesh(displacement);
    std::vector<double> box2 = mesh->getBoundingBox();
    displacement = std::vector<double>(mesh->getDim(),-1.0);
    double volume = 1.0;
    for (int i=0; i<mesh->getDim(); i++)
        volume *= box1[2*i+1]-box1[2*i+0];
    if ( volume > 0.0 )
        utils->passes("non-zero bounding box");
    else
        utils->failure("non-zero bounding box");
    bool passes = true;
    for (size_t i=0; i<box1.size(); i++) {
        if ( fabs(box2[i]-box1[i]-1.0) > 1e-12 )
            passes = false;
    }
    if ( passes )
        utils->passes("scalar displacement test");
    else
        utils->failure("scalar displacement test");
    // Test displacement vector
    #ifdef USE_AMP_VECTORS
        // Get the volume of each element
        size_t numElements = mesh->numLocalElements(mesh->getGeomType());
        std::vector<double> orig_vol(numElements,0.0);
        AMP::Mesh::MeshIterator  cur_elem = mesh->getIterator(mesh->getGeomType(),0);
        for (size_t i=0; i<numElements; i++) {
            orig_vol[i] = cur_elem->volume();
            ++cur_elem;
        }
        // Get the position of the nodes
        AMP::LinearAlgebra::Vector::shared_ptr  posVec1 = mesh->getPositionVector ( "pos_before" );
        // Displace the mesh
        AMP::LinearAlgebra::Vector::shared_ptr  dispVec = posVec1->cloneVector ( "displ" );
        dispVec->copyVector(posVec1);
        dispVec->scale(1e-3);
        mesh->displaceMesh ( dispVec );
        // Get the new positions
        AMP::LinearAlgebra::Vector::shared_ptr  posVec2 = mesh->getPositionVector ( "pos_after" );
        AMP::LinearAlgebra::Vector::shared_ptr diff = dispVec->cloneVector ( "diff" );
        diff->subtract ( posVec2 , posVec1 );
        if ( diff->equals ( dispVec ) )
            utils->passes ( "displacement successfully applied" );
        else
            utils->failure ( "displacement failed" );
        // Get the new volumes
        bool volume_passed = true;
        cur_elem = mesh->getIterator(mesh->getGeomType(),0);
        double vol_ratio = 1.0;
        for (int i=0; i<mesh->getDim(); i++)
            vol_ratio *= 1.001;
        for (size_t i=0; i<numElements; i++) {
            double ratio = cur_elem->volume() / orig_vol[i];
            if ( !AMP::Utilities::approx_equal(ratio,vol_ratio,1e-9) )      // The new volume should be (1+10^-3)^dim the original volume
                volume_passed = false;
            ++cur_elem;
        }
        if ( volume_passed )
            utils->passes ( "displacement changed volumes" );
        else
            utils->failure ( "displacement changed volumes" );
    #endif
}


// Test getting parent elements for each mesh element
void getParents( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) 
{
    bool pass = true;
    int gcw = mesh->getMaxGhostWidth();
    for (int type1=0; type1<=(int)mesh->getGeomType(); type1++) {
        AMP::Mesh::MeshIterator it = mesh->getIterator((AMP::Mesh::GeomType)type1,gcw);
        for (size_t k=0; k<it.size(); k++) {
            for (int type2=0; type2<type1; type2++) {
                std::vector<AMP::Mesh::MeshElement> elements = it->getElements((AMP::Mesh::GeomType)type2);
                for (size_t i=0; i<elements.size(); i++) {
                    if ( !elements[i].globalID().is_local() )
                        continue;
                    std::vector<AMP::Mesh::MeshElement> parents = mesh->getElementParents(elements[i],(AMP::Mesh::GeomType)type1);
                    // Check that the current parent was found (find all parents)
                    bool found = false;
                    for (size_t j=0; j<parents.size(); j++) {
                        if ( parents[j]==*it )
                            found = true;
                    }
                    if ( !found )
                        pass = false;
                    // Check that all parents do have the current element as a child (no extra parents found)
                    for (size_t j=0; j<parents.size(); j++) {
                        std::vector<AMP::Mesh::MeshElement> children = parents[j].getElements((AMP::Mesh::GeomType)type2);
                        found = false;
                        for (size_t m=0; m<children.size(); m++) {
                            if ( children[m]==elements[i] )
                                found = true;
                        }
                        if ( !found )
                            pass = false;
                    }
                }
            }
            ++it;
        }
    }
    if ( pass )
        utils->passes("getParents passed");
    else
        utils->failure("getParents passed");    
}



/*

class  VerifyElementForNode
{
public:
    static const char * get_test_name () { return "verify beginElementForNode"; }

    static  bool element_has_node ( AMP::Mesh::MeshAdapter::Element e , AMP::Mesh::MeshAdapter::Node n ) {
        for ( size_t i = 0 ; i != e.numNodes(); i++ )
        {
            if ( n.globalID() == e.getNodeID ( i ) )
            {
              return true;
            }
        }
        return false;
    }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
          AMP::Mesh::MeshAdapter::OwnedNodeIterator curOwned = mesh->beginOwnedNode();
          bool passedTest = true;
          while ( curOwned != mesh->endOwnedNode() )
          {
            AMP::Mesh::MeshAdapter::NodeElementIterator  curElem = mesh->beginElementForNode ( *curOwned );
            if ( !element_has_node ( *curElem , *curOwned ) )
            {
              passedTest = false;
            }
            curOwned++;
          }
          if ( passedTest )
            utils->passes ( "All elements found are correct" );
          else
            utils->failure ( "Found an incorrect element" );
        }
    };

    struct ElementHelper
    {
        typedef AMP::Mesh::MeshAdapter::ElementIterator  Iterator;
        static Iterator  begin ( AMP::Mesh::MeshAdapter::shared_ptr &p ) { return p->beginElement(); }
        static Iterator  end ( AMP::Mesh::MeshAdapter::shared_ptr &p ) { return p->endElement(); }
    };
    struct NodeHelper
    {
        typedef AMP::Mesh::MeshAdapter::NodeIterator  Iterator;
        static Iterator  begin ( AMP::Mesh::MeshAdapter::shared_ptr &p ) { return p->beginNode(); }
        static Iterator  end ( AMP::Mesh::MeshAdapter::shared_ptr &p ) { return p->endNode(); }
};


template <typename HELPER>
class  VerifyProcAndIsOwnedInterface
{
public:
    static const char * get_test_name () { return "Verify isOwned Interface"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
        std::map<size_t , std::vector<size_t> > Ids;
        typename HELPER::Iterator cur = HELPER::begin ( mesh );
        bool totalTest = true;
        while ( cur != HELPER::end ( mesh ) )
        {
            Ids[cur->procID()].push_back ( cur->globalID() );
            if ( (int)cur->procID() == utils->rank() )
            {
              if ( !cur->isOwned() )
                totalTest = false;
            }
            else
            {
              if ( cur->isOwned() )
                totalTest = false;
            }
            cur++;
        }
        if ( totalTest )
            utils->passes ( "isOwned works" );
        else
            utils->failure ( "isOwned fails" );

        totalTest = true;
        for (int i=0; i!=utils->size(); i++ )
        {
            size_t res = Ids[i].size();
            res = globalComm.bcast(res,i);
            if ( res > 0 )
            {
              std::vector<size_t>  incoming;
              if ( i == utils->rank() )
              {
                globalComm.bcast(&(Ids[i][0]),(int)res,i);
              }
              else
              {
                incoming.resize ( res );
                globalComm.bcast(&(incoming[0]),(int)res,i);
                for (size_t j=0; j!=Ids[i].size(); i++ )
                {
                  bool found = false;
                  for ( size_t k = 0 ; k != incoming.size() ; k++ )
                  {
                    if ( Ids[i][j] == incoming[k] )
                    {
                      found = true;
                      break;
                    }
                  }
                  if ( !found )
                    totalTest = false;
                }
              }
            }
        }
        if ( totalTest )
            utils->passes ( "procID is set correctly" );
        else
            utils->failure ( "procID is not set correctly" );
    }
};



AMP::Mesh::MeshAdapter::shared_ptr globalMeshForMeshVectorFactory = AMP::Mesh::MeshAdapter::shared_ptr();
template <int SIZE, bool NODAL, bool RUNTIME>
class  MeshVectorFactory
{
public:
    typedef  AMP::LinearAlgebra::Vector               vector;

    static  AMP::LinearAlgebra::Variable::shared_ptr  getVariable ()
    {
        if ( NODAL ) {
            return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,SIZE>( "test vector" ) );
        } else {
            if ( RUNTIME )
                return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::Mesh::RunTimeIntegrationPointVariable ( "test vector" , SIZE ) );
            else
                return AMP::LinearAlgebra::Variable::shared_ptr ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable,SIZE>( "test vector" ) );
        }
    }

    static  AMP::LinearAlgebra::Vector::shared_ptr getVector()
    {
        if ( globalMeshForMeshVectorFactory.get()==NULL )
            AMP_ERROR("mesh must be set before this can be called");
        return globalMeshForMeshVectorFactory->createVector ( getVariable() );
    }

    static  AMP::Mesh::DOFMap::shared_ptr getDOFMap()
    {
        if ( globalMeshForMeshVectorFactory.get()==NULL )
            AMP_ERROR("mesh must be set before this can be called");
        return globalMeshForMeshVectorFactory->getDOFMap ( getVariable() );
    }

};


// Test bug when subsetting a multivector for a variable
class Bug_758
{
public:
        static const char * get_test_name () { return "Bug_758"; }

        static  bool verify_memory_address ( AMP::LinearAlgebra::Vector::shared_ptr  multi , AMP::LinearAlgebra::Vector::shared_ptr first , AMP::LinearAlgebra::Vector::shared_ptr second )
        {
          multi->setRandomValues ();
          AMP::LinearAlgebra::Vector::iterator  cur_multi , cur_sub;
          cur_multi = multi->begin();
          cur_sub = first->begin();
          bool retVal = true;
          size_t i = 0;
          while ( cur_multi != multi->end() )
          {
            if ( cur_sub == first->end() )
            {
              cur_sub = second->begin();
            }
            if ( cur_sub == second->end() )
            {
              retVal = false;
              break;
            }
            if ( *cur_sub != *cur_multi )
            {
              retVal = false;
              break;
            }
            i++;
            cur_multi++;
            cur_sub++;
          }
          if ( i != multi->getLocalSize() )
          {
            retVal = false;
          }
          if ( cur_sub != second->end() )
          {
            retVal = false;
          }
          return retVal;
        }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
          AMP::LinearAlgebra::Variable::shared_ptr var1 ( new AMP::Mesh::NodalScalarVariable ( "number_1" ) );
          AMP::LinearAlgebra::Variable::shared_ptr var2 ( new AMP::Mesh::Nodal3VectorVariable ( "number_2" ) );
          AMP::LinearAlgebra::Variable::shared_ptr multivar ( new AMP::LinearAlgebra::MultiVariable ( "multi" ) );
          multivar->castTo<AMP::LinearAlgebra::MultiVariable>().add ( var1 );
          multivar->castTo<AMP::LinearAlgebra::MultiVariable>().add ( var2 );

          AMP::LinearAlgebra::Vector::shared_ptr  p1 = mesh->createVector ( multivar );
          AMP::LinearAlgebra::Vector::shared_ptr  p2 = p1->subsetVectorForVariable ( var1 );
          AMP::LinearAlgebra::Vector::shared_ptr  p3 = p1->subsetVectorForVariable ( var2 );
          if ( *(p2->getVariable()) == *var1 )
            utils->passes ( "First variable is correct" );
          else
            utils->failure ( "First variable is incorrect" );
          if ( *(p3->getVariable()) == *var2 )
            utils->passes ( "Second variable is correct" );
          else
            utils->failure ( "Second variable is incorrect" );
          if ( verify_memory_address ( p1 , p2 , p3 ) )
            utils->passes ( "Subset vector pulls out the correct vectors 1" );
          else
            utils->failure ( "Subset vector fails to pull out the correct vectors 1" );
          AMP::LinearAlgebra::Vector::shared_ptr  p4 = AMP::LinearAlgebra::PetscVector::view ( p1 );
          AMP::LinearAlgebra::Vector::shared_ptr  p5 = p4->subsetVectorForVariable ( var1 );
          AMP::LinearAlgebra::Vector::shared_ptr  p6 = p4->subsetVectorForVariable ( var2 );
          if ( *(p5->getVariable()) == *var1 )
            utils->passes ( "First variable is correct" );
          else
            utils->failure ( "First variable is incorrect" );
          if ( *(p6->getVariable()) == *var2 )
            utils->passes ( "Second variable is correct" );
          else
            utils->failure ( "Second variable is incorrect" );
          if ( verify_memory_address ( p4 , p5 , p6 ) )
            utils->passes ( "Subset vector pulls out the correct vectors 2" );
          else
            utils->failure ( "Subset vector fails to pull out the correct vectors 2" );

          Vec vec = p4->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
          AMP::LinearAlgebra::Vector::shared_ptr p7 ( reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *> ( vec->data ) , AMP::LinearAlgebra::ExternalVectorDeleter() );
          AMP::LinearAlgebra::Vector::shared_ptr p8 = p7->subsetVectorForVariable ( var1 );
          AMP::LinearAlgebra::Vector::shared_ptr p9 = p7->subsetVectorForVariable ( var2 );
          if ( *(p8->getVariable()) == *var1 )
            utils->passes ( "First variable is correct" );
          else
            utils->failure ( "First variable is incorrect" );
          if ( *(p9->getVariable()) == *var2 )
            utils->passes ( "Second variable is correct" );
          else
            utils->failure ( "Second variable is incorrect" );
          if ( verify_memory_address ( p7 , p8 , p9 ) )
            utils->passes ( "Subset vector pulls out the correct vectors 3" );
          else
            utils->failure ( "Subset vector fails to pull out the correct vectors 3" );

        }
};



template <int SIZE>
class Bug_761
{
public:
        static const char * get_test_name () { return "Verify fix to Bug 761"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        AMP::LinearAlgebra::Variable::shared_ptr  var ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::IntegrationPointVariable , SIZE> ( "t" ) );
        AMP::Mesh::DOFMap::shared_ptr  dofmap = mesh->getDOFMap ( var );
        AMP::LinearAlgebra::Vector::shared_ptr  vector = mesh->createVector ( var );

        vector->setToScalar ( 123.4 );
        bool testPass = true;
        AMP::Mesh::MeshAdapter::ElementIterator curElem = mesh->beginElement();
        while ( curElem != mesh->endElement() )
        {
            for ( size_t i = 0 ; i != SIZE ; i++ )
            {
              if ( vector->getValueByGlobalID ( dofmap->getGlobalID ( curElem->globalID() , i ) ) != 123.4 )
              {
                testPass = false;
              }
            }
            curElem++;
        }
        if ( testPass )
            utils->passes ( "Every value of a integration point vector found" );
        else
            utils->failure ( "Failed to find an integration point vector" );
    }
};


class Bug_623
{
public:
        static const char * get_test_name () { return "Verify fix to Bug 623"; }

    static  void run_test ( AMP::UnitTest *utils ) {
          AMP::Mesh::MeshAdapter::shared_ptr  mesh1 ( new AMP::Mesh::LibMeshAdapter () );
          AMP::Mesh::MeshAdapter::shared_ptr  mesh2 ( new AMP::Mesh::LibMeshAdapter () );
          AMP::LinearAlgebra::Variable::shared_ptr     var ( new AMP::Mesh::NodalScalarVariable ( "temp" ) );

          // Create new meshes
          mesh1->generateCube ( 5 );
          mesh2->generateCube ( 5 );

          // On mesh1, create matrix then vector
          mesh1->createMatrix ( var );
          mesh1->createVector ( var );

          // On mesh2, create vector then matrix
          mesh2->createVector ( var );
          mesh2->createMatrix ( var );
          utils->passes ( "Bug 623" );
        }
};


class  VerifyNodeElemMapIteratorTest
{
public:
        static const char * get_test_name () { return "Verify Node<->Element iterator"; }

        static bool  verify_node ( AMP::Mesh::MeshAdapter::Node n , AMP::Mesh::MeshAdapter::shared_ptr mesh )
        {
          std::set<size_t>  elems_from_node , elems_from_mesh;
          AMP::Mesh::MeshAdapter::NodeElementIterator  cur_elem = mesh->beginElementForNode ( n );
          while ( cur_elem != mesh->endElementForNode ( n ) )
          {
            elems_from_node.insert ( cur_elem->globalID() );
            cur_elem++;
          }
          AMP::Mesh::MeshAdapter::ElementIterator  cur_elem2 = mesh->beginElement ();
          while ( cur_elem2 != mesh->endElement() )
          {
            for ( size_t i = 0 ; i != cur_elem2->numNodes() ; i++ )
              if ( cur_elem2->getNodeID(i) == n.globalID() )
                elems_from_mesh.insert ( cur_elem2->globalID() );
            cur_elem2++;
          }

          if ( elems_from_node.size() == 0 ) return false;
          if ( std::equal ( elems_from_node.begin() , elems_from_node.end() , elems_from_mesh.begin() ) )
              return true;
          return false;
        }


    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {

          int i = 0;
          int SKIP = (utils->size() * mesh->numTotalNodes()) / 20;

          AMP::Mesh::MeshAdapter::NodeIterator  cur_node = mesh->beginNode();
          while ( cur_node != mesh->endNode() )
          {
            if ( i == SKIP )
            {
              if ( !verify_node ( *cur_node , mesh ) )
              {
                utils->failure ("Verify Node<->Element iterator");
                return;
              }
              i = 0;
            }
            else
              i++;
            cur_node++;
          }

          utils->passes ("Verify Node<->Element iterator");
        }
};


class VerifyBoundaryIteratorTest
{
public:
        static const char * get_test_name () { return "Verify BoundaryIterator"; }

        static int  isBoundaryElement ( AMP::Mesh::MeshAdapter::Element e )
        {
          for ( size_t i = 0 ; i != e.numSides() ; i++ )
            if ( !e.hasNeighbor ( i ) ) return true;
          return false;
        }

        static bool isNodeOnBoundary ( AMP::Mesh::MeshAdapter::Node  n ,  AMP::Mesh::MeshAdapter::shared_ptr  mesh )
        {
          AMP::Mesh::MeshAdapter::NodeElementIterator  cur_elem = mesh->beginElementForNode ( n );
          while ( cur_elem != mesh->endElementForNode ( n ) )
          {
            if ( !isBoundaryElement ( *cur_elem ) )
              return false;
            cur_elem++;
          }
          return true;
        }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {

          std::set<short int>::const_iterator  cur_bid = mesh->getBoundaryIds().begin();
          while ( cur_bid != mesh->getBoundaryIds().end() )
          {
            AMP::Mesh::MeshAdapter::BoundaryNodeIterator  cur_b_node = mesh->beginBoundary ( *cur_bid );
            while ( cur_b_node != mesh->endBoundary ( *cur_bid ) )
            {
              if ( !isNodeOnBoundary ( *cur_b_node , mesh ) )
              {
                utils->failure ("Verify Node<->Element iterator");
                return;
              }
              cur_b_node++;
            }
            cur_bid++;
          }
          utils->passes ("Verify Node<->Element iterator");
        }
};


*/

#endif

