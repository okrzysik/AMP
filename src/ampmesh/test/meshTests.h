#ifndef included_MeshTests
#define included_MeshTests

#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/MeshIterator.h"


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
        if ( id.is_local )
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
        if ( id.is_local ) {
            for (int i=0; i<=(int)type; i++) {
                if ( i!=0 && i!=(int)type )
                    continue;  // getElements is unfinished for types other than verticies and elements
                AMP::Mesh::GeomType type2 = (AMP::Mesh::GeomType) i;
                std::vector<AMP::Mesh::MeshElement> pieces = element.getElements(type2);
            }
            std::vector< AMP::Mesh::MeshElement::shared_ptr > neighbors = element.getNeighbors();
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
            // Check that the mesh has ghost cells if gcw>0 and N_proc>1
            if ( gcw>0 && mesh->getComm().getSize()>1 ) {
                if ( type==0 || type==mesh->getGeomType() ) {
                    if ( mesh->numGhostElements(type,gcw) > 0 )
                        ut->passes("Mesh has ghost cells");
                    else 
                        ut->failure("Mesh has ghost cells");
                } else {
                    ut->expected_failure("ghost cell test with intermediate element type is not supported yet");
                }
            }
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


// Test the number of elements in the mesh
void MeshCountTest( AMP::UnitTest *ut, boost::shared_ptr<AMP::Mesh::Mesh> mesh )
{
    AMP::AMP_MPI comm = mesh->getComm();
    for (int i=0; i<=(int)mesh->getGeomType(); i++) {
        AMP::Mesh::GeomType type = (AMP::Mesh::GeomType) i;
        if ( type!=AMP::Mesh::Vertex && type!=mesh->getGeomType() ) {
            ut->expected_failure("Not finished testing all element type in MeshCountTest");
            continue;
        }
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
        if ( comm.getSize() > 1 && type!=0 ) {
            if ( mesh->numGhostElements(type,1)!=0 && (mesh->numGhostElements(type,1)+N_local)!=N_global )
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
    size_t meshID = mesh->meshID();
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
}



// This tests loops over all boundary ids
void VerifyBoundaryNodeIterator( AMP::UnitTest *utils, AMP::Mesh::Mesh::shared_ptr mesh ) {
    const std::vector<int> bids = mesh->getIDSets();
    for (size_t i=0; i<bids.size(); i++) {
        int bid = bids[i];
        for (int gcw=0; gcw<=1; gcw++) {
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
            if ( node_ids.size()==0 )
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


/*


class  DisplaceNodes
{
public:
    static const char * get_test_name () { return "displace nodes"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        // Get the volume of each element
        std::vector<double> orig_vol(mesh->numLocalElements(),0.0);
        AMP::Mesh::MeshAdapter::ElementIterator  cur_elem = mesh->beginElement();
        for (size_t i=0; i<mesh->numLocalElements(); i++) {
            orig_vol[i] = (*cur_elem).volume();
            cur_elem++;
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
        bool volume_changed = true;
        cur_elem = mesh->beginElement();
        for (size_t i=0; i<mesh->numLocalElements(); i++) {
            double vol_ratio = (*cur_elem).volume() / orig_vol[i];
            if ( !AMP::Utilities::approx_equal(vol_ratio,1.003003001,1e-9) )        // The new volume should be (1+10^-3)^dim the original volume
                volume_changed = false;
            cur_elem++;
        }
        if ( volume_changed )
            utils->passes ( "displacement changed volumes" );
        else
            utils->failure ( "displacement changed volumes" );
    }
};


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

