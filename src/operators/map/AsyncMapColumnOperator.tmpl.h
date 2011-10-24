
#include "operators/map/AsyncMapOperator.h"
#include "operators/map/AsyncMapOperatorParameters.h"
#include "operators/map/Map3to1to3.h"

namespace AMP {
namespace Operator {


  template <typename MAP_TYPE>
  boost::shared_ptr<AsyncMapColumnOperator>  AsyncMapColumnOperator::build ( AMP::Mesh::MeshManager::shared_ptr meshManager , boost::shared_ptr<Database> db )
  {
    typedef boost::shared_ptr < AsyncMapOperator >             AsyncOp_ptr;
    typedef boost::shared_ptr < AsyncMapOperatorParameters >   AsyncOpParams_ptr;
    Map3to1to3::spMap::type       sharedMap ( new std::multimap<double,double> );
    sharedMap.isComputed = false;

    std::vector < AsyncOp_ptr >           asyncToDo;
    std::vector < AsyncOpParams_ptr >     asyncParams;

    boost::shared_ptr<AsyncMapColumnOperatorParameters>  newParams ( new AsyncMapColumnOperatorParameters ( db ) );
    boost::shared_ptr<AsyncMapColumnOperator>  newMapColumn ( new AsyncMapColumnOperator ( newParams ) );

    for ( size_t i = 0 ; i != meshManager->getNumMaps() ; i++ )
    {
      if ( !MAP_TYPE::validMapType ( meshManager->getMapType ( i ) ) )
      {
        continue;
      }

      AMP::Mesh::MeshManager::Adapter::shared_ptr  mesh = meshManager->getMesh ( meshManager->getMapMeshName ( i ) );
      AsyncOpParams_ptr  mapParams ( new typename MAP_TYPE::Parameters ( meshManager->getMapDB ( i ) ) );

      mapParams->d_MapComm           = meshManager->getMapComm ( i );
      mapParams->d_SetBegin          = mesh->beginOwnedBoundary ( meshManager->getMapBoundaryId ( i ) );
      mapParams->d_SetEnd            = mesh->endOwnedBoundary ( meshManager->getMapBoundaryId ( i ) );
      mapParams->d_IsMaster          = meshManager->getMapDominance ( i ) == AMP::Mesh::MeshManager::Master;
      mapParams->d_MeshAdapter       = mesh;
      mapParams->d_ConstructionPhase = 0;

      size_t  tagOffset = meshManager->getMapTagOffset ( i );
      mapParams->d_ToMasterCommTag = MAP_TYPE::CommTagBase + tagOffset;
      mapParams->d_ToSlaveCommTag  = MAP_TYPE::CommTagBase + tagOffset + 1;
      if ( meshManager->getMapConstructionParam ( i ) == AMP::Mesh::MeshManager::Synchronous )
      {
        mapParams->d_AsynchronousConstructionParam  = 0;
      }
      else
      {
        if ( mapParams->d_IsMaster )
        {
          mapParams->d_AsynchronousConstructionParam  = 1;
        }
        else
        {
          mapParams->d_AsynchronousConstructionParam  = 2;
        }
      }

      AsyncOp_ptr  mapOp ( new MAP_TYPE ( mapParams ) );
      boost::shared_ptr<Map3to1to3>  curM313;
      if ( curM313 = boost::dynamic_pointer_cast<Map3to1to3> ( mapOp ) )
      {
        curM313->getMapData() = sharedMap;
      }

      newMapColumn->append ( mapOp );

      if ( meshManager->getMapConstructionParam ( i ) == AMP::Mesh::MeshManager::Asynchronous )
      {
        asyncToDo.push_back ( mapOp );
        asyncParams.push_back ( mapParams );
      }
    }


    while ( asyncToDo.size() )
    {
      std::vector<AsyncOp_ptr>::iterator        curOp;
      std::vector<AsyncOpParams_ptr>::iterator  curParams;

      curOp = asyncToDo.begin();
      curParams = asyncParams.begin();

      while ( curOp != asyncToDo.end() )
      {
        if ( (*curOp)->continueAsynchronousConstruction ( *curParams ) )
        {
          curOp = asyncToDo.erase ( curOp );
          curParams = asyncParams.erase ( curParams );
        }
        else
        {
          curOp++;
          curParams++;
        }
      }
    }
    return newMapColumn;
  }

}
}

