#include "ampmesh/MeshAdapter.h"
#include "ampmesh/LibMeshAdapter.h"
#include "ampmesh/MeshManager.h"
#include "../contact_search/ContactManagerTmpl.h"
#include "test_MeshGenerators.h"

namespace AMP {
namespace unit_test {


template <typename GENERATOR=LibMeshCubeGenerator<5> >
class InstantiateCoarseSearch
{
public:
    static const char * get_test_name () 
    { 
        return "Instantiate a CoarseSearch object"; 
    }

    static  void run_test ( UnitTest *utils )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        typedef AMP::Mesh::CoarseSearch<AMP::Mesh::MeshManager::Adapter , AMP::Mesh::MeshManager::Adapter::NodeIterator , AMP::Mesh::MeshManager::Adapter::ElementIterator>  CS;
        CS::Parameters::shared_ptr params ( 
                  new CS::Parameters ( mesh->beginNode() ,
                                       mesh->endNode () , 
                                       mesh->beginElement () ,
                                       mesh->endElement () ) );
        CS coarseSearch ( params );
        utils->passes ( "CoarseSearch instantiated" );
    }
};


template <typename GENERATOR=LibMeshCubeGenerator<5> >
class VerifyContactSearchSerialPenaltyInterface
{
public:
    static const char * get_test_name ()
    {
        return "Verify Contact Search Serial Penalty Interaface";
    }

    static  void run_test ( UnitTest *utils )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        AMP::Mesh::SerialPenaltyContactSearch::Side1Set              s1Set;
        AMP::Mesh::SerialPenaltyContactSearch::Side2Set              s2Set;
        AMP::Mesh::SerialPenaltyContactSearch::Side1ToSide2MultiMap  s1Map;
        AMP::Mesh::SerialPenaltyContactSearch::Side2ToSide1MultiMap  s2Map;

        AMP::Mesh::SerialPenaltyContactSearch::computeContactSets ( mesh ,
                                             mesh ,
                                             mesh->beginBoundary ( 1 ) ,
                                             mesh->endBoundary ( 1 ) ,
                                             mesh->beginSideBoundary ( 1 ) ,
                                             mesh->endSideBoundary ( 1 ) ,
                                             mesh->beginSideBoundary ( 1 ) ,
                                             mesh->endSideBoundary ( 1 ) ,
                                             s1Set ,
                                             s2Set ,
                                             s1Map ,
                                             s2Map );
          
        AMP::Mesh::SerialPenaltyContactSearch::Side2ToSide1MultiMap::iterator  curPair = s2Map.begin();
        AMP::Mesh::SerialPenaltyContactSearch::Side2Iterator                   curFacet = curPair->first;
        std::set<size_t>                                            nodesFound;

        bool failed = false;
        while ( curPair != s2Map.end() ) {
            if ( curPair->first != curFacet ) {
                if ( nodesFound.size() != 4 ) {
                    utils->failure ( "Found the wrong number of nodes" );
                    failed = true;
                }
                for ( size_t i = 0 ; i != 4 ; i++ ) {
                    if ( nodesFound.find ( curFacet->getNodeID ( i ) ) == nodesFound.end() ) {
                        utils->failure ( "Did not find a corner node!" );
                        failed = true;
                    }
                }
                nodesFound.clear();
                curFacet = curPair->first;
                if ( failed ) {
                    break;
                }
            }
            nodesFound.insert ( curPair->second->globalID() );
            curPair++;
        }
        if ( !failed )
            utils->passes ( "SerialPenaltyContactSearch finds nodes on corners" );
    }
};


template <typename GENERATOR=LibMeshCubeGenerator<5> >
class CallContactSearchSerialPenaltyInterface
{
public:
    static const char * get_test_name ()
    {
        return "Call Contact Search Serial Penalty Interaface";
    }

    static  void run_test ( UnitTest *utils )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        AMP::Mesh::SerialPenaltyContactSearch::Side1Set              s1Set;
        AMP::Mesh::SerialPenaltyContactSearch::Side2Set              s2Set;
        AMP::Mesh::SerialPenaltyContactSearch::Side1ToSide2MultiMap  s1Map;
        AMP::Mesh::SerialPenaltyContactSearch::Side2ToSide1MultiMap  s2Map;

        AMP::Mesh::SerialPenaltyContactSearch::computeContactSets ( mesh ,
                                             mesh ,
                                             mesh->beginBoundary ( 1 ) ,
                                             mesh->endBoundary ( 1 ) ,
                                             mesh->beginSideBoundary ( 1 ) ,
                                             mesh->endSideBoundary ( 1 ) ,
                                             mesh->beginSideBoundary ( 1 ) ,
                                             mesh->endSideBoundary ( 1 ) ,
                                             s1Set ,
                                             s2Set ,
                                             s1Map ,
                                             s2Map );
        utils->passes ( "computeContactSets executes on a trivial set" );
    }
};


template <typename GENERATOR=LibMeshCubeGenerator<5> >
class SimpleCoarseSearch
{
public:
    static const char * get_test_name () 
    { 
        return "Simple coarse search test of same mesh with same iterators"; 
    }

    static  void run_test ( UnitTest *utils )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        typedef AMP::Mesh::CoarseSearch<AMP::Mesh::MeshManager::Adapter , 
                                    AMP::Mesh::MeshManager::Adapter::ElementIterator , 
                                    AMP::Mesh::MeshManager::Adapter::ElementIterator>  CS;
        CS::Parameters::shared_ptr params ( 
                  new CS::Parameters ( mesh->beginElement() ,
                                       mesh->endElement () , 
                                       mesh->beginElement () ,
                                       mesh->endElement () ) );
        params->d_mesh1 = mesh;
        params->d_mesh2 = mesh;
        CS coarseSearch ( params );
        typename CS::FirstToSecondMap  m1ToM2;
        typename CS::SecondToFirstMap  m2ToM1;

        coarseSearch.findNeighbors ( m1ToM2 , m2ToM1 );
        if ( m1ToM2.size() > 0 )
            utils->passes ( "Maps are nontrivial" );
        else
            utils->failure ( "Maps are trivial" );

        if ( m1ToM2.size() == m2ToM1.size() )
            utils->passes ( "Maps are the same size" );
        else
            utils->failure ( "Maps are different sizes" );

        typename CS::FirstToSecondMap::iterator  curItem = m1ToM2.begin();
        AMP::Mesh::MeshManager::Adapter::ElementIterator  curElem = curItem->first;
        bool  foundMyself = false;;
        while ( curItem != m1ToM2.end() ) {
            if ( curItem->first != curElem ) {
                if ( !foundMyself ) {
                    utils->failure ( "Didn't find myself in coarse search" );
                    break;
                }
                foundMyself = false;
                curElem = curItem->first;
            }
            size_t id1 = curItem->second->globalID();
            size_t id2 = curElem->globalID();
            if ( id1 == id2 )
                foundMyself = true;
            typename CS::SecondToFirstMap::iterator  curMapped = m2ToM1.lower_bound ( curItem->second );
            typename CS::SecondToFirstMap::iterator  curMappedEnd = m2ToM1.upper_bound ( curItem->second );
            bool foundCurrentMap = false;
            while ( curMapped != curMappedEnd ) {
                if ( curMapped->second->globalID() == curElem->globalID() ) {
                    foundCurrentMap = true;
                    break;
                }
                curMapped++;
            }
            if ( !foundCurrentMap ) {
                utils->failure ( "Maps are not inverse" );
                break;
            }
            curItem++;
        }
        utils->passes ( "Iterated through maps" );
    }
};


template <typename GENERATOR=LibMeshCubeGenerator<5> >
class FindNodesInElementsCoarseSearch
{
public:
    static const char * get_test_name () 
    { 
        return "Attempt to find all nodes belonging to an element"; 
    }

    static  void run_test ( UnitTest *utils )
    {
        GENERATOR meshGenerator;
        AMP::Mesh::MeshAdapter::shared_ptr mesh = meshGenerator.getMesh();
        typedef AMP::Mesh::CoarseSearch<AMP::Mesh::MeshManager::Adapter , 
                                    AMP::Mesh::MeshManager::Adapter::NodeIterator , 
                                    AMP::Mesh::MeshManager::Adapter::ElementIterator>  CS;
        CS::Parameters::shared_ptr params ( 
                  new CS::Parameters ( mesh->beginNode() ,
                                       mesh->endNode () , 
                                       mesh->beginElement () ,
                                       mesh->endElement () ) );
        params->d_mesh1 = mesh;
        params->d_mesh2 = mesh;
        CS coarseSearch ( params );
        typename CS::FirstToSecondMap  m1ToM2;
        typename CS::SecondToFirstMap  m2ToM1;

        coarseSearch.findNeighbors ( m1ToM2 , m2ToM1 );
        if ( m1ToM2.size() > 0 )
            utils->passes ( "Maps are nontrivial" );
        else
            utils->failure ( "Maps are trivial" );

        typename CS::SecondToFirstMap::iterator  curItem = m2ToM1.begin();
        AMP::Mesh::MeshManager::Adapter::ElementIterator  curElem = curItem->first;
        std::set<size_t>  nodes;
        for ( size_t i = 0 ; i != curElem->numNodes() ; i++ )
            nodes.insert ( curElem->getNodeID ( i ) );
        size_t numElems = 1;
        while ( curItem != m2ToM1.end() ) {
            if ( curElem != curItem->first ) {
                curElem = curItem->first;
                if ( nodes.size() > 0 ) {
                    utils->failure ( "Didn't find nodes on an element" );
                    break;
                }
                for ( size_t i = 0 ; i != curElem->numNodes() ; i++ )
                    nodes.insert ( curElem->getNodeID ( i ) );
                numElems++;
            }
            std::set<size_t>::iterator loc = nodes.find ( curItem->second->globalID() );
            if ( loc != nodes.end() ) {
                nodes.erase ( loc );
            }
            curItem++;
        }
        utils->passes ( "Iterated through maps" );

        curElem = mesh->beginElement();
        while ( curElem != mesh->endElement() ) {
            curElem++;
            numElems--;
        }
        if ( numElems == 0 )
            utils->passes ( "Found all nodes on all elements" );
        else
            utils->passes ( "Missed some elements" );
    }
};


}
}

