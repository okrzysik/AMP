#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/AMP_MPI.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Writer.h"
#include "utils/Utilities.h"
#include "utils/PIO.h"
#include "ampmesh/Mesh.h"
#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "operators/map/ScalarN2GZAxisMap.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/ColumnOperator.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include "discretization/createLibmeshElements.h"

#define __INIT_FN1__(x, y, z) ( x+y+z )
#define __INIT_FN2__(x, y, z) ( 2*x+y+z )
#define __INIT_FN3__(x, y, z) ( 4*x+y+z )

int runTest(std::string exeName, AMP::UnitTest *ut)
{
  // Input and output file names
  std::string input_file =  exeName;
  std::string log_file = "output_" + exeName;

  AMP::PIO::logAllNodes(log_file);
  boost::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  // Construct a mesh manager which reads in the fuel mesh
  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  boost::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  boost::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr manager(AMP::Mesh::Mesh::buildMesh(mgrParams) );
  AMP::pout << "Finished loading meshes" <<  std::endl;

  int nodalGhostWidth = 0;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr phiDofMap      = AMP::Discretization::simpleDOFManager::create(manager, AMP::Mesh::Vertex, nodalGhostWidth,      1,    split);
  AMP::Discretization::DOFManager::shared_ptr eectDofMap     = AMP::Discretization::simpleDOFManager::create(manager, AMP::Mesh::Vertex, nodalGhostWidth,      5,    split);

  //Construct Variable
  AMP::LinearAlgebra::Variable::shared_ptr potentialVariable(new AMP::LinearAlgebra::Variable("Potential"));
  AMP::LinearAlgebra::Vector::shared_ptr   potentialMapVec          = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, true);
  AMP::LinearAlgebra::Vector::shared_ptr   potentialSolVec          = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, true);
  AMP::LinearAlgebra::Vector::shared_ptr   potentialResVec          = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, true);

  //Construct Variable
  boost::shared_ptr<AMP::LinearAlgebra::Variable> batteryVariables (new AMP::LinearAlgebra::Variable("Battery"));
  AMP::LinearAlgebra::Vector::shared_ptr BatterySolVec              = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , true);
  AMP::LinearAlgebra::Vector::shared_ptr BatteryResidualVec         = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , true);
  AMP::LinearAlgebra::Vector::shared_ptr BatteryMapVec              = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , true);

  AMP::LinearAlgebra::Vector::shared_ptr ElectrodeSolVec            = BatterySolVec->select( AMP::LinearAlgebra::VS_Stride( 3, 5) , "V4" );
  AMP::LinearAlgebra::Vector::shared_ptr ElectrodeMapVec            = BatteryMapVec->select( AMP::LinearAlgebra::VS_Stride( 3, 5) , "V4" );
//---------------------------------------------------

  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
  siloWriter->registerMesh( manager );
  siloWriter->setDecomposition(1);
  siloWriter->registerVector( potentialMapVec , manager, AMP::Mesh::Vertex, "potentialVec" );
  siloWriter->registerVector( BatterySolVec , manager, AMP::Mesh::Vertex, "batteryVec" );
  siloWriter->writeFile( log_file , 0);

//---------------------------------------------------

  boost::shared_ptr<AMP::LinearAlgebra::MultiVector> multiSolutionVec;
  multiSolutionVec = AMP::LinearAlgebra::MultiVector::create( "MultiSolutionVec", globalComm );
  multiSolutionVec->addVector(ElectrodeSolVec);
  multiSolutionVec->addVector(potentialSolVec);

  boost::shared_ptr<AMP::LinearAlgebra::MultiVector> multiResVec;
  multiResVec = AMP::LinearAlgebra::MultiVector::create( "MultiRHSVec", globalComm );
  multiResVec->addVector(BatteryResidualVec);
  multiResVec->addVector(potentialResVec);

   // Make new vectors
  boost::shared_ptr<AMP::LinearAlgebra::MultiVector> multiSolutionMapVec;
  multiSolutionMapVec = AMP::LinearAlgebra::MultiVector::create( "MultiSolutionMapVec", globalComm );
  multiSolutionMapVec->addVector(ElectrodeMapVec);
  multiSolutionMapVec->addVector(potentialMapVec);

  ////=-------------------------------------------------------------
  // make map operator
  boost::shared_ptr<AMP::Database> mapOperatorDatabase = input_db->getDatabase("PotentialMaps");
  boost::shared_ptr<AMP::Operator::Map3to1to3Parameters> mapOperatorParameters1(new AMP::Operator::Map3to1to3Parameters(mapOperatorDatabase));
  mapOperatorParameters1->d_Mesh1= manager->Subset("CellSandwich_2_1");
  mapOperatorParameters1->d_Mesh2= manager->Subset("AnodeCC_1_1");
  mapOperatorParameters1->d_BoundaryID1= 5;
  mapOperatorParameters1->d_BoundaryID2= 1;
  mapOperatorParameters1->d_MapComm = globalComm;
//  boost::shared_ptr<AMP::Operator::StridedZAxisMap> mapOperator1(new AMP::Operator::StridedZAxisMap(mapOperatorParameters1));

//  mapOperator1->setVector ( multiSolutionMapVec );
//  mapOperator1->apply(multiRhsVec, multiSolutionVec, multiResVec);

  boost::shared_ptr<AMP::Operator::Map3to1to3Parameters> mapOperatorParameters2(new AMP::Operator::Map3to1to3Parameters(mapOperatorDatabase));
  mapOperatorParameters2->d_Mesh1= manager->Subset("CellSandwich_2_1");
  mapOperatorParameters2->d_Mesh2= manager->Subset("CathodeCC_3_1");
  mapOperatorParameters2->d_BoundaryID1= 3;
  mapOperatorParameters2->d_BoundaryID2= 2;
  mapOperatorParameters2->d_MapComm = globalComm;
//  boost::shared_ptr<AMP::Operator::StridedZAxisMap> mapOperator2(new AMP::Operator::StridedZAxisMap(mapOperatorParameters2));

//  mapOperator2->setVector ( multiSolutionMapVec );
//  mapOperator2->apply(multiRhsVec, multiSolutionVec, multiResVec);

  boost::shared_ptr<AMP::Operator::ColumnOperator> potentialMapsColumn;
  boost::shared_ptr<AMP::Operator::ColumnOperatorParameters>  mapColParams ( new AMP::Operator::ColumnOperatorParameters ( input_db ) );
  potentialMapsColumn.reset( new AMP::Operator::ColumnOperator ( mapColParams ) );
//  potentialMapsColumn->append( mapOperator1 );
//  potentialMapsColumn->append( mapOperator2 );

  AMP::Mesh::MeshIterator node  = anodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
  AMP::Mesh::MeshIterator end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN1__(pt[0], pt[1], pt[2]);

        potentialSolVec->setValueByGlobalID(bndGlobalIds[0], val);
  }//end for node

  node  = cathodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN3__(pt[0], pt[1], pt[2]);

        potentialSolVec->setValueByGlobalID(bndGlobalIds[0], val);
  }//end for node

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 3, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN2__(pt[0], pt[1], pt[2]);

        ElectrodeSolVec->setValueByGlobalID(bndGlobalIds[3], val);
  }//end for node

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 5, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN2__(pt[0], pt[1], pt[2]);

        ElectrodeSolVec->setValueByGlobalID(bndGlobalIds[3], val);
  }//end for node
  multiSolutionVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  potentialMapsColumn->apply(multiRhsVec, multiSolutionVec, multiResVec);

  node  = anodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val1 = __INIT_FN2__(pt[0], pt[1], pt[2]);
        double val2 = potentialMapVec->getValueByGlobalID(bndGlobalIds[0]);

        if ( !AMP::Utilities::approx_equal(val1,val2) )
        {
            ut->passes(" DTK Map Operator anodeCCtest ");
        } else {
            ut->failure(" DTK Map Operator test ");
        }
  }//end for node

  node  = cathodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val1 = __INIT_FN2__(pt[0], pt[1], pt[2]);
        double val2 = potentialMapVec->getValueByGlobalID(bndGlobalIds[0]);

        if ( !AMP::Utilities::approx_equal(val1,val2) )
        {
            ut->passes(" DTK Map Operator cathodeCCtest ");
        } else {
            ut->failure(" DTK Map Operator test ");
        }
  }//end for node

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 3, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val1 = __INIT_FN1__(pt[0], pt[1], pt[2]);
        double val2 = BatteryMapVec->getValueByGlobalID(bndGlobalIds[3]);

        if ( !AMP::Utilities::approx_equal(val1,val2) )
        {
            ut->passes(" DTK Map Operator CellSandwich test ");
        } else {
            ut->failure(" DTK Map Operator test ");
        }
  }//end for node

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 5, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val1 = __INIT_FN3__(pt[0], pt[1], pt[2]);
        double val2 = BatteryMapVec->getValueByGlobalID(bndGlobalIds[3]);

        if ( !AMP::Utilities::approx_equal(val1,val2) )
        {
            ut->passes(" DTK Map Operator CellSandwich test ");
        } else {
            ut->failure(" DTK Map Operator test ");
        }
  }//end for node
}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    AMP::AMP_MPI globalComm = AMP::AMP_MPI(AMP_COMM_WORLD);
    //int  numNodes = globalComm.getSize();
    runTest ( "input_testDTKMapOperator-1" , &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
