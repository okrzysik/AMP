#ifndef included_MeshTests
#define included_MeshTests

#include "utils/AMP_MPI.h"
#include "utils/UnitTest.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "ampmesh/MeshIterator.h"


// This test checks that we can iterate over different elements
// Note: in this context, elements refer to the highest dimension object
class  ElementIteratorTest
{
public:
    static const char * get_test_name () { return "verify element iterators"; }

    static  void run_test ( AMP::UnitTest *utils, boost::shared_ptr<AMP::Mesh::Mesh> mesh ) {
        // Get the number of elements
        AMP::Mesh::GeomType type = mesh->getGeomType();
        int number_of_local_elements = (int) mesh->numLocalElements(type);
        if ( number_of_local_elements > 0 )
            utils->passes ( "non trivial mesh generated" );
        else
            utils->failure ( "trivial mesh generated" );
        // Get the element iterator
        AMP::Mesh::MeshIterator  cur_elem = mesh->getIterator(type,0);
        AMP::Mesh::MeshIterator  end_elem = cur_elem.end();
        std::set<size_t>  ids;
        while ( cur_elem != end_elem )
        {
            ids.insert ( cur_elem->globalID() );
            number_of_local_elements--;
            cur_elem++;
        }
        if ( number_of_local_elements == 0 )
            utils->passes ( "regular iterator count" );
        else
            utils->failure ( "regular iterator count" );
        if ( ids.size() == mesh->numLocalElements(type) )
            utils->passes ( "regular iterator uniqueness" );
        else
            utils->failure ( "regular iterator uniqueness" );

/*        number_of_local_elements = mesh->numLocalElements ();
        AMP::Mesh::MeshAdapter::ConstElementIterator  cur_const_elem = mesh->beginElement();
        ids.clear();
        while (!( cur_const_elem == mesh->endElement() ))
        {
            ids.insert ( cur_const_elem->globalID() );
            number_of_local_elements--;
            cur_const_elem++;
        }
        if ( number_of_local_elements == 0 )
            utils->passes ( "const iterator" );
        else
            utils->failure ( "const iterator" );
        if ( ids.size() == mesh->numLocalElements () )
            utils->passes ( "const iterator uniqueness" );
        else
            utils->failure ( "const iterator uniqueness" );
    */
    }
};

/*
class  ElementTest
{
public:
    static const char * get_test_name () { return "check elements"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        AMP::Mesh::MeshAdapter::ElementIterator  cur_elem = mesh->beginElement();
        bool volume_pass = true;
        while ( cur_elem != mesh->endElement() ) {
            if ( (*cur_elem).volume() <= 0.0 )
                volume_pass = false;
            cur_elem++;
        }
        if ( volume_pass )
            utils->passes ( "elements have volumes > 0" );
        else
            utils->failure ( "elements have volumes > 0" );
    }
};


class  VerifyOwnedBoundaryNodeIterator
{
public:
    static const char * get_test_name () { return "verfify BoundaryNodeIterator interface"; }

    static  void verify_boundary ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh, short int bid ) {
        AMP::Mesh::MeshAdapter::OwnedBoundaryNodeIterator curNode = mesh->beginOwnedBoundary ( bid );
        std::set<size_t>  node_ids;
        while ( curNode != mesh->endOwnedBoundary ( bid ) )
        {
            node_ids.insert ( curNode->globalID() );
            curNode++;
        }

        bool testPassed = true;
        size_t  numFound = 0;
        AMP::Mesh::MeshAdapter::OwnedNodeIterator  curMNode = mesh->beginOwnedNode();
        while ( curMNode != mesh->endOwnedNode() )
        {
            if ( mesh->isOnBoundary ( *curMNode , bid ) )
            {
              numFound++;
              if ( node_ids.find ( curMNode->globalID() ) == node_ids.end() )
              {
                testPassed = false;
                break;
              }
            }
            curMNode++;
        }
        if ( testPassed )
            utils->passes ( "Found all boundary nodes" );
        else
            utils->failure ( "Node on boundary not in iterator" );
        if ( numFound == node_ids.size() )
            utils->passes ( "All boundary nodes accounted for" );
        else
            utils->failure ( "Node in iterator not on boundary" );
    }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        const std::set<short int> bids = mesh->getBoundaryIds ();
        std::set<short int>::const_iterator curBid = bids.begin();
        while ( curBid != bids.end() )
        {
            verify_boundary ( utils, mesh, *curBid );
            curBid++;
        }
    }
};


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


class  VerifyOwnedNodeIterator
{
public:
        static const char * get_test_name () { return "verify owned node iterator"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        AMP::Mesh::MeshAdapter::OwnedNodeIterator   curONode = mesh->beginOwnedNode ();
        std::set<size_t>  ids1,ids2;
        while ( curONode != mesh->endOwnedNode () )
        {
            ids1.insert ( curONode->globalID() );
            curONode++;
        }
        AMP::Mesh::MeshAdapter::NodeIterator  curNode = mesh->beginNode();
        while ( curNode != mesh->endNode() )
        {
            if ( curNode->isOwned() )
            {
              ids2.insert ( curNode->globalID() );
            }
            curNode++;
        }
        if ( ids1.size() == ids2.size() )
        {
            if ( std::equal ( ids1.begin() , ids1.end() , ids2.begin() ) )
              utils->passes ( "All owned nodes found" );
            else
              utils->failure ( "Owned node sets are different sizes" );
        }
        else
        {
            if ( ids1.size() > ids2.size() )
              utils->passes ( "OwnedNodeIterator found too many nodes" );
            else
              utils->passes ( "OwnedNodeIterator did not find enough nodes" );
        }
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


class  VerifyBoundaryNodeIterator
{
public:
    static const char * get_test_name () { return "verfify BoundaryNodeIterator interface"; }

    static  void verify_boundary ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh, short int bid ) {
        AMP::Mesh::MeshAdapter::BoundaryNodeIterator curNode = mesh->beginBoundary ( bid );
        std::set<size_t>  node_ids;
        while ( curNode != mesh->endBoundary ( bid ) )
        {
            size_t  gid = curNode->globalID();
            node_ids.insert ( gid );
            curNode++;
        }

        bool testPassed = true;
        AMP::Mesh::MeshAdapter::NodeIterator  curMNode = mesh->beginNode();
        while ( curMNode != mesh->endNode() )
          {
            if ( mesh->isOnBoundary ( *curMNode , bid ) )
            {
              if ( node_ids.find ( curMNode->globalID() ) == node_ids.end() )
              {
                testPassed = false;
                break;
              }
              else
              {
                node_ids.erase ( curMNode->globalID() );
              }
            }
            curMNode++;
        }
        if ( testPassed )
            utils->passes ( "Found all boundary nodes" );
        else
            utils->failure ( "Node on boundary not in iterator" );

        /** It would be nice to be able to do this test, but libmesh insists on preventing me.
        size_t k = node_ids.size();
        if ( k == 0 )
            utils->passes ( "All boundary nodes accounted for" );
        else
        {
            std::stringstream msg;
            msg << "Node " << (*(node_ids.begin())) << " in iterator not on boundary!";
            utils->failure ( msg.str().c_str() );
        }
            */ /*
    }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
        const std::set<short int> bids = mesh->getBoundaryIds ();
        std::set<short int>::const_iterator curBid = bids.begin();
        while ( curBid != bids.end() )
        {
            verify_boundary ( utils, mesh, *curBid );
            curBid++;
        }
    }
};


class  NodeIteratorTest
{
public:
        static const char * get_test_name () { return "verify node iterators"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {
          int number_of_local_nodes = mesh->numLocalNodes ();
          AMP::Mesh::MeshAdapter::NodeIterator  cur_node = mesh->beginNode();
          std::set<size_t>  ids;
          while ( cur_node != mesh->endNode() )
          {
            ids.insert ( cur_node->globalID() );
            number_of_local_nodes--;
            cur_node++;
          }
          if ( number_of_local_nodes == 0 )
            utils->passes ( "regular iterator count" );
          else
            utils->failure ( "regular iterator count" );
          if ( ids.size() == mesh->numLocalNodes () )
            utils->passes ( "regular iterator uniqueness" );
          else
            utils->failure ( "regular iterator uniqueness" );


          number_of_local_nodes = mesh->numLocalNodes ();
          AMP::Mesh::MeshAdapter::ConstNodeIterator  cur_const_node = mesh->beginNode();
          ids.clear();
          while (!( cur_const_node == mesh->endNode() ))
          {
            ids.insert ( cur_const_node->globalID() );
            number_of_local_nodes--;
            cur_const_node++;
          }
          if ( number_of_local_nodes == 0 )
            utils->passes ( "const iterator" );
          else
            utils->failure ( "const iterator" );
          if ( ids.size() == mesh->numLocalNodes () )
            utils->passes ( "const iterator uniqueness" );
          else
            utils->failure ( "const iterator uniqueness" );
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


template <int DOF_PER_NODE>
class VerifyGetVectorTest
{
public:
    static const char * get_test_name () { return "Verify vector interface in MeshAdapter"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {

        AMP::LinearAlgebra::Variable::shared_ptr variable ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,DOF_PER_NODE> ( "test vector" ) );
        AMP::LinearAlgebra::Vector::shared_ptr vectora = mesh->createVector ( variable );  // Generates new vector
        AMP::LinearAlgebra::Vector::shared_ptr vectorb = mesh->createVector ( variable );  // Gets from the cached copy

        size_t  num_dofs = mesh->numTotalNodes() * DOF_PER_NODE;
        if ( vectora->getGlobalSize() == num_dofs )
            utils->passes ( "global vector size" );
        else
            utils->failure ( "global vector size" );

        vectora->setRandomValues ();
        double t1 = vectora->L2Norm();
        vectora->abs ( vectora );
        if ( fabs ( vectora->L2Norm() - t1 ) < 0.0000001 )
            utils->passes ( "non-trivial random vector" );
        else
            utils->failure ( "non-trivial random vector" );

        vectorb->setToScalar ( 3. );
        vectora->multiply ( vectora , vectorb );
        if ( fabs ( vectora->L2Norm() - 3.*t1 ) < 0.00000001 )
            utils->passes ( "trivial usage" );
        else
            utils->failure ( "trivial usage" );

        // Verify some math...
        globalMeshForMeshVectorFactory = mesh;
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,false,false> > ( utils );
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,false,true > > ( utils );
        test_managed_vectors_loop< MeshVectorFactory<DOF_PER_NODE,true, false> > ( utils );
        test_parallel_vectors_loop<MeshVectorFactory<DOF_PER_NODE,true, false> > ( utils );
        globalMeshForMeshVectorFactory = AMP::Mesh::MeshAdapter::shared_ptr();
    }
};


template <int DOF_PER_NODE>
class VerifyGetMatrixTrivialTest
{
public:
    static const char * get_test_name () { return "Verify matrix interface trivially in MeshAdapter"; }

    static  void run_test ( AMP::UnitTest *utils, AMP::Mesh::MeshAdapter::shared_ptr mesh ) {

          AMP::LinearAlgebra::Variable::shared_ptr variable ( new AMP::LinearAlgebra::VectorVariable<AMP::Mesh::NodalVariable,DOF_PER_NODE> ( "test vector" ) );
          AMP::LinearAlgebra::Matrix::shared_ptr matrixa = mesh->createMatrix ( variable );  // Generates new mtrix
          AMP::LinearAlgebra::Vector::shared_ptr vectorb = mesh->createVector ( variable );  // Gets vector from the cached copy
          AMP::LinearAlgebra::Vector::shared_ptr vectorc = mesh->createVector ( variable );  // Gets vector from the cached copy

          vectorb->setRandomValues ();
          matrixa->makeConsistent ();
          matrixa->mult ( vectorb , vectorc );
          if ( vectorc->L1Norm() < 0.00000001 )
            utils->passes ( "obtained 0 matrix from mesh" );
          else
            utils->failure ( "did not obtain 0 matrix from mesh" );

          AMP::LinearAlgebra::Matrix::shared_ptr matrixb = mesh->createMatrix ( variable );   // Need to get another matrix to
                                                                           // store data due to Epetra insert/replace
                                                                           // idiom.  Matrixa is fixed with no entires.
          vectorc->setToScalar ( 1. );
          matrixb->makeConsistent ();
          matrixb->setDiagonal ( vectorc );
          matrixb->mult ( vectorb , vectorc );
          vectorb->subtract ( vectorb , vectorc );

          if ( vectorb->L1Norm() < 0.0000001 )
            utils->passes ( "created identity matrix from mesh" );
          else
            utils->failure ( "created identity matrix from mesh" );
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

