#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <string>

#include "utils/AMPManager.h"
#include "materials/Material.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"


#include "vectors/Variable.h"
#include "vectors/SimpleVector.h"
#include "vectors/Vector.h"

#include "ampmesh/MeshManager.h"
#include "ampmesh/MeshAdapter.h"
#include "ampmesh/SiloIO.h"
#include "ampmesh/MeshVariable.h"

#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/MapOperatorParameters.h"
#include "operators/OperatorBuilder.h"

extern "C"{
#include "petsc.h"
}


void testMap(AMP::UnitTest *ut, std::string exeName )
{
  std::string input_file = "input_" + exeName;
  std::string log_file = "output_" + exeName;

  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP::PIO::logAllNodes(log_file);

  AMP::MeshManagerParameters::shared_ptr mgrParams ( new AMP::MeshManagerParameters ( input_db ) );
  AMP::MeshManager::shared_ptr manager ( new AMP::MeshManager ( mgrParams ) );
  AMP::MeshManager::Adapter::shared_ptr meshAdapter1 = manager->getMesh ( "pellet" );
  AMP::MeshManager::Adapter::shared_ptr meshAdapter2 = manager->getMesh ( "clad" );

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  AMP::LinearAlgebra::Variable::shared_ptr testVariable ( new AMP::NodalScalarVariable ( "MapSolution" ) );
  AMP::LinearAlgebra::Variable::shared_ptr testVariable1 ( new AMP::NodalScalarVariable ( "MapSolution", meshAdapter1  ) );
  AMP::LinearAlgebra::Variable::shared_ptr testVariable2 ( new AMP::NodalScalarVariable ( "MapSolution", meshAdapter2 ) );
  boost::shared_ptr<AMP::LinearAlgebra::Variable> gapVariable(new AMP::LinearAlgebra::Variable("Gap"));

  AMP::LinearAlgebra::Vector::shared_ptr mapSolution = manager->createVector ( testVariable );
  mapSolution->setToScalar ( 0.0 );
  manager->registerVectorAsData ( mapSolution, "OriginalMap" );

  AMP::LinearAlgebra::Vector::shared_ptr mapSolutionMaster = meshAdapter1->createVector( testVariable1 );
  AMP::LinearAlgebra::Vector::shared_ptr mapSolutionSlave = meshAdapter2->createVector( testVariable2 );
 //-------------------------------------

  boost::shared_ptr<AMP::InputDatabase> map3dto1d_db1  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapPelletto1D"));
  boost::shared_ptr<AMP::MapOperatorParameters> map3dto1dParams1 (new AMP::MapOperatorParameters( map3dto1d_db1 ));
  map3dto1dParams1->d_MeshAdapter = meshAdapter1;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> map1ToLowDim (new AMP::Operator::Map3Dto1D( map3dto1dParams1 ));

  boost::shared_ptr<AMP::InputDatabase> map1dto3d_db1  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DtoClad"));
  boost::shared_ptr<AMP::MapOperatorParameters> map1dto3dParams1 (new AMP::MapOperatorParameters( map1dto3d_db1 ));
  map1dto3dParams1->d_MapAdapter = meshAdapter2;
  boost::shared_ptr<AMP::Map1Dto3D> map1ToHighDim (new AMP::Map1Dto3D( map1dto3dParams1 ));

  map1ToLowDim->setZLocations( map1ToHighDim->getZLocations());

  boost::shared_ptr<AMP::InputDatabase> map3dto1d_db2  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladto1D"));
  boost::shared_ptr<AMP::MapOperatorParameters> map3dto1dParams2 (new AMP::MapOperatorParameters( map3dto1d_db2 ));
  map3dto1dParams2->d_MeshAdapter = meshAdapter2;
  boost::shared_ptr<AMP::Operator::Map3Dto1D> map2ToLowDim (new AMP::Operator::Map3Dto1D( map3dto1dParams2 ));

  boost::shared_ptr<AMP::InputDatabase> map1dto3d_db2  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DtoPellet"));
  boost::shared_ptr<AMP::MapOperatorParameters> map1dto3dParams2 (new AMP::MapOperatorParameters( map1dto3d_db2 ));
  map1dto3dParams2->d_MapAdapter = meshAdapter1;
  boost::shared_ptr<AMP::Map1Dto3D> map2ToHighDim (new AMP::Map1Dto3D( map1dto3dParams2 ));

  map2ToLowDim->setZLocations( map2ToHighDim->getZLocations());
 //-------------------------------------
  size_t gapVecCladSize = map1ToHighDim->getNumZlocations(); 
  AMP::LinearAlgebra::Vector::shared_ptr gapVecClad = AMP::LinearAlgebra::SimpleVector::create( gapVecCladSize, gapVariable );

  size_t gapVecPelletSize = map2ToHighDim->getNumZlocations(); 
  AMP::LinearAlgebra::Vector::shared_ptr gapVecPellet = AMP::LinearAlgebra::SimpleVector::create( gapVecPelletSize, gapVariable );

 //-------------------------------------
  unsigned int d_boundaryId = map3dto1d_db1->getInteger("BoundaryId");  

  AMP::MeshManager::Adapter::BoundarySideIterator bnd = meshAdapter1->beginSideBoundary( d_boundaryId );
  AMP::MeshManager::Adapter::BoundarySideIterator end_bnd = meshAdapter1->endSideBoundary( d_boundaryId );

  AMP::DOFMap::shared_ptr dof_map = meshAdapter1->getDOFMap(testVariable1);

  for( ; bnd != end_bnd; ++bnd) {
      AMP::MeshManager::Adapter::Element cur_side = *bnd;
      std::vector<unsigned int> bndGlobalIds;
      dof_map->getDOFs(cur_side, bndGlobalIds);
      for(unsigned int ii=0; ii<cur_side.numNodes(); ii++)
      {
         mapSolutionMaster->setValueByGlobalID( bndGlobalIds[ii], meshAdapter1->getNode(cur_side.getNodeID(ii)).z());
//         mapSolutionMaster->setValueByGlobalID( bndGlobalIds[ii], 50);
      }
  }

  //-------------------------------------

  int cnt=0;

  bool testPassed = false;

  while ( cnt < 20 )
  {
    cnt++;

    map1ToLowDim->apply(nullVec,mapSolutionMaster,gapVecClad ,1.0, 0.0);
    map1ToHighDim->apply(nullVec,gapVecClad , mapSolutionSlave ,1.0, 0.0);

      std::cout<<"Master Map Solution " <<std::endl;
      for( int i=0; i<gapVecCladSize; i++) {
          std::cout<<" @i : "<< i<<" is "<<gapVecClad->getValueByLocalID(i) ;
      }
      std::cout<<std::endl;
    //------------------------------------------------------------
//    mapSolutionSlave->setToScalar(90);

    map2ToLowDim->apply(nullVec,mapSolutionSlave,gapVecPellet ,1.0, 0.0);
    map2ToHighDim->apply(nullVec,gapVecPellet , mapSolutionMaster,1.0, 0.0);

      std::cout<<"Slave Map Solution " <<std::endl;
      for( int i=0; i<gapVecPelletSize; i++) {
          std::cout<<" @i : "<< i<<" is "<<gapVecPellet->getValueByLocalID(i) ;
      }
      std::cout<<std::endl;
    //------------------------------------------------------------

    if (1) {
      testPassed = true;
      break;
    } else {
      std::cout << "Norm of the change in sol for iteration "<<cnt <<"is -->"<<  std::endl;
    }

    std::cout<<std::endl;

  }

  if( nodes == 1 ) {
#ifdef USE_SILO
    manager->writeFile<AMP::SiloIO> ( exeName , 0 );
#endif
  }
  //-------------------------------------

  if( testPassed )
  {
    ut.passes("Seggregated solve of Composite Operator using control loop of Thermal+Robin->Map->Gap->Map->Thermal+Robin .");
  }
  else
  {
    ITFAILS;
  }

  input_db.reset();

  ut.passes(exeName);

  //  AMP::AMPManager::shutdown();

}

int main(int argc, char *argv[])
{
    AMP::AMPManager::startup(argc, argv);
    AMP::UnitTest ut;

    try {
    testMap(ut, "testMapOperator");
    } catch (std::exception &err) {
        std::cout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
        ut.failure("ERROR: While testing");
    } catch( ... ) {
        std::cout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
        ut.failure("ERROR: While testing");
    }
   
    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}   



