#ifndef included_MeshTests
#define included_MeshTests

#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/SubsetMesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/MeshIterator.h"

#ifdef USE_AMP_VECTORS
    #include "vectors/Vector.h"
#endif


// This test checks a single mesh element iterator
// ut           Unit test class to report the results
// iterator     local iterator over elements
// N_local      number of local elements for the iterator
// N_ghost      number of ghost elements for the iterator
void ElementIteratorTest( AMP::UnitTest *ut, AMP::Mesh::MeshIterator iterator, 
    const size_t N_local, const size_t N_ghost, const AMP::Mesh::GeomType type )
{
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

    // Run element tests
    bool id_pass = true;
    bool type_pass = true;
    bool volume_pass = true;
    bool coord_pass = true;
    bool elements_pass = true;
    bool neighbor_pass = true;
    cur_it = iterator.begin();
    while ( cur_it != end_it ) {
        AMP::Mesh::MeshElement element = *cur_it;
        AMP::Mesh::MeshElementID id = element.globalID();
        if ( id != cur_it->globalID() )
            id_pass = false;
        if ( element.elementType() != type )
            type_pass = false;
        if ( type==AMP::Mesh::Vertex ) {
            std::vector<double> coord = element.coord();
            if ( coord.size()==0 )
                coord_pass = false;
        } else {
            if ( element.volume() <= 0.0 )
                volume_pass = false;
        }
        if ( id.is_local() ) {
            for (int i=0; i<=(int)type; i++) {
                AMP::Mesh::GeomType type2 = (AMP::Mesh::GeomType) i;
                std::vector<AMP::Mesh::MeshElement> pieces = element.getElements(type2);
                if ( pieces.empty() )
                    elements_pass = false;
            }
            std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = element.getNeighbors();
            if ( neighbors.empty() )
                neighbor_pass = false;
        }
        ++cur_it;   // Pre-increment is faster than post-increment
    }
    if ( id_pass && type_pass && volume_pass && coord_pass ) {
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
        if ( !elements_pass )
            ut->failure( "elements failed getElements test" );
        if ( !neighbor_pass )
            ut->failure( "elements failed getNeighbors test" );
    }
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
                ElementIteratorTest( ut, iterator, N_local, N_ghost, type );
            // Add tests with gcw != 0
            // Add const iterator tests
        }
    }    
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
        size_t N_sum = comm.sumReduce(N_local);
        if ( N_global > 0 )
            ut->passes("Non-trival mesh created");
        else
            ut->failure("Non-trival mesh created");
        if ( N_sum == N_global )
            ut->passes("Sum of local mesh counts matches global count");
        else
            ut->failure("Sum of local mesh counts matches global count");
        if ( mesh->numGhostElements(type,0) == 0 )
            ut->passes("gcw=0 has no ghost elements");
        else
            ut->failure("gcw=0 has no ghost elements");
    }
}


// Test some basic Mesh properties
void MeshBasicTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    // test that we can get the mesh ID
    AMP::Mesh::MeshID meshID = mesh->meshID();
    if ( meshID>0 && meshID<1000 ) 
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
}


// This tests loops over all boundary ids
void VerifyBoundaryIDNodeIterator( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) {
    const std::vector<int> bids = mesh->getIDSets();
    for (size_t i=0; i<bids.size(); i++) {
        int bid = bids[i];
        for (int gcw=0; gcw<=0; gcw++) {
            // Get the iterator over the current boundary id
            AMP::Mesh::MeshIterator curNode = mesh->getIDsetIterator( AMP::Mesh::Vertex, bid, gcw );
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
void VerifyBoundaryIterator( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) {
    for (int gcw=0; gcw<=0; gcw++) {
        for (int type2=0; type2<=(int)mesh->getGeomType(); type2++) {
            AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) type2;
            // Get the iterator over the current boundary id
            AMP::Mesh::MeshIterator iterator = mesh->getSurfaceIterator( type, gcw );
            size_t global_size = mesh->getComm().sumReduce(iterator.size());
            bool passes = global_size>0;
            if ( boost::dynamic_pointer_cast<AMP::Mesh::SubsetMesh>(mesh).get()==NULL )
                passes = passes && global_size<mesh->numGlobalElements(type);
            if ( passes )
                utils->passes("Non-trivial surface iterator created");
            else
                utils->failure("Non-trivial surface iterator created");
        }
    }
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
    // Test the default values
    if ( id0.meshID()!=static_cast<size_t>(-1) || id0.is_local() || id0.type()!=0 || id0.owner_rank()!=0 || id0.local_id()!=static_cast<unsigned int>(-1) )
        utils->failure("MeshElementID test defaults");
    if ( num_failed0 == utils->NumFailLocal() )
        utils->passes("MeshElementID tests");
    else
        utils->failure("MeshElementID tests");
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

