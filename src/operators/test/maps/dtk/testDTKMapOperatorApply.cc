#include <ampmesh/Mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <operators/ColumnOperator.h>
#include <operators/map/dtk/DTKMapOperator.h>
#include <operators/map/dtk/MultiDofDTKMapOperator.h>
#include <utils/AMPManager.h>
#include <utils/AMP_MPI.h>
#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/PIO.h>
#include <utils/UnitTest.h>
#include <utils/Utilities.h>
#include <utils/Writer.h>
#include <vectors/MultiVector.h>
#include <vectors/Variable.h>
#include <vectors/VectorBuilder.h>
#include <vectors/VectorSelector.h>

#

void thermalTest(AMP::UnitTest *ut, std::string input_file)
{
  std::string log_file = "log_DTKMapOperatorApply" ;
  std::string out_file = "out_DTKMapOperatorApply";

  AMP::PIO::logOnlyNodeZero(log_file);

  AMP::shared_ptr<AMP::InputDatabase> input_db(new AMP::InputDatabase("input_db"));
  AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

  //PROFILE_START("SetupDriver");
  AMP::InputManager::getManager()->parseInputFile(input_file, input_db);
  input_db->printClassData(AMP::plog);

  AMP_INSIST(input_db->keyExists("Mesh"), "Key ''Mesh'' is missing!");
  AMP::shared_ptr<AMP::Database>  mesh_db = input_db->getDatabase("Mesh");
  AMP::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(new AMP::Mesh::MeshParameters(mesh_db));
  mgrParams->setComm(AMP::AMP_MPI(AMP_COMM_WORLD));
  AMP::Mesh::Mesh::shared_ptr manager(AMP::Mesh::Mesh::buildMesh(mgrParams) );
  AMP::pout << "Finished loading meshes" <<  std::endl;

  int DOFsPerNode = 1;
  int nodalGhostWidth = 1;
  int DOFsPerElement = 8;
  int gaussPointGhostWidth = 1;
  bool split = true;
  AMP::Discretization::DOFManager::shared_ptr nodalDofMap      = AMP::Discretization::simpleDOFManager::create(manager, AMP::Mesh::Vertex, nodalGhostWidth,      DOFsPerNode,    split);

  AMP::LinearAlgebra::Variable::shared_ptr thermalVariable(new AMP::LinearAlgebra::Variable("Temperature"));

  AMP::LinearAlgebra::Vector::shared_ptr  SolutionVec           = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
  AMP::LinearAlgebra::Vector::shared_ptr  RightHandSideVec      = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
  AMP::LinearAlgebra::Vector::shared_ptr  ResidualVec           = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable );
  
  AMP::LinearAlgebra::Vector::shared_ptr  thermalMapVec         = AMP::LinearAlgebra::createVector( nodalDofMap, thermalVariable , true);

  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
  siloWriter->registerMesh( manager );
  siloWriter->registerVector( SolutionVec , manager, AMP::Mesh::Vertex, "SolutionVec" );
  siloWriter->registerVector( thermalMapVec, manager, AMP::Mesh::Vertex, "MapVec" );

  RightHandSideVec->setToScalar(0.0);

  std::vector<AMP::Mesh::MeshID> meshIDs = manager->getBaseMeshIDs();
 
  siloWriter->writeFile( out_file , 0);

  AMP::shared_ptr<AMP::Operator::ColumnOperator> thermalMapsColumn;

  AMP::Mesh::Mesh::shared_ptr cellSandwichMesh = manager->Subset("CellSandwich");
  AMP::Mesh::Mesh::shared_ptr ccMesh           = manager->Subset("CellCurrentCollectors");

  SolutionVec->setToScalar(298.);           
  thermalMapVec->setToScalar(-1);           
  RightHandSideVec->setToScalar(0);           
      

  AMP::pout<<"----------------------------\n";
  AMP::pout<<"     CREATE MAP OPERATOR    \n";
  AMP::pout<<"----------------------------\n";
  AMP::shared_ptr<AMP::Database> nullDatabase;
  // INTERFACE WITH ANODE
  AMP::pout<<"interface anodeCC cellSandwich\n";
  AMP::shared_ptr<AMP::Operator::MultiDofDTKMapOperatorParameters> anodeCCCellSandwichMapOperatorParams(new AMP::Operator::MultiDofDTKMapOperatorParameters(nullDatabase));
  anodeCCCellSandwichMapOperatorParams->d_globalComm    = AMP_COMM_WORLD;
  anodeCCCellSandwichMapOperatorParams->d_Mesh1 = cellSandwichMesh;
  anodeCCCellSandwichMapOperatorParams->d_BoundaryID1 = 1;
  anodeCCCellSandwichMapOperatorParams->d_Variable1 = thermalVariable ->getName();
  anodeCCCellSandwichMapOperatorParams->d_StrideOffset1 = 0;
  anodeCCCellSandwichMapOperatorParams->d_StrideLength1 = 1;
  anodeCCCellSandwichMapOperatorParams->d_Mesh2 = ccMesh;
  anodeCCCellSandwichMapOperatorParams->d_BoundaryID2 = 3;
  anodeCCCellSandwichMapOperatorParams->d_Variable2 = thermalVariable ->getName();
  anodeCCCellSandwichMapOperatorParams->d_StrideOffset2 = 0;
  anodeCCCellSandwichMapOperatorParams->d_StrideLength2 = 1;
  anodeCCCellSandwichMapOperatorParams->d_SourceVector = SolutionVec;
  anodeCCCellSandwichMapOperatorParams->d_TargetVector = thermalMapVec;
  AMP::shared_ptr<AMP::Operator::Operator> anodeCCCellSandwichMapOperator(new AMP::Operator::MultiDofDTKMapOperator(anodeCCCellSandwichMapOperatorParams));

  // INTERFACE WITH CATHODE
  AMP::pout<<"interface cellSandwich cathodeCC\n";
  AMP::shared_ptr<AMP::Operator::MultiDofDTKMapOperatorParameters> cellSandwichCathodeCCMapOperatorParams(new AMP::Operator::MultiDofDTKMapOperatorParameters(nullDatabase));
  cellSandwichCathodeCCMapOperatorParams->d_globalComm    = AMP_COMM_WORLD;
  cellSandwichCathodeCCMapOperatorParams->d_Mesh1 = cellSandwichMesh;
  cellSandwichCathodeCCMapOperatorParams->d_BoundaryID1 = 2;
  cellSandwichCathodeCCMapOperatorParams->d_Variable1 = thermalVariable ->getName();
  cellSandwichCathodeCCMapOperatorParams->d_StrideOffset1 = 0;
  cellSandwichCathodeCCMapOperatorParams->d_StrideLength1 = 1;
  cellSandwichCathodeCCMapOperatorParams->d_Mesh2 = ccMesh;
  cellSandwichCathodeCCMapOperatorParams->d_BoundaryID2 = 4;
  cellSandwichCathodeCCMapOperatorParams->d_Variable2 = thermalVariable ->getName();
  cellSandwichCathodeCCMapOperatorParams->d_StrideOffset2 = 0;
  cellSandwichCathodeCCMapOperatorParams->d_StrideLength2 = 1;
  cellSandwichCathodeCCMapOperatorParams->d_SourceVector = SolutionVec;
  cellSandwichCathodeCCMapOperatorParams->d_TargetVector = thermalMapVec;
  AMP::shared_ptr<AMP::Operator::Operator> cellSandwichCathodeCCMapOperator(new AMP::Operator::MultiDofDTKMapOperator(cellSandwichCathodeCCMapOperatorParams));


  AMP::shared_ptr<AMP::Operator::ColumnOperatorParameters>  mapColParams ( new AMP::Operator::ColumnOperatorParameters ( input_db ) );
  thermalMapsColumn.reset( new AMP::Operator::ColumnOperator ( mapColParams ) );
  thermalMapsColumn->append( anodeCCCellSandwichMapOperator  );
  thermalMapsColumn->append( cellSandwichCathodeCCMapOperator);

  AMP::LinearAlgebra::Vector::shared_ptr nullVec;

  thermalMapsColumn->apply(SolutionVec, ResidualVec);
  AMP::pout<< " L2Norm of Map Vec " << std::setprecision(17)   << thermalMapVec->L2Norm() << std::endl;

  siloWriter->writeFile( out_file , 1);
}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  std::string inputFile = "input_testDTKMapOperatorApply";

  try {
    if(argc>1) 
      inputFile = argv[1];
    thermalTest(&ut, inputFile );
  } catch (std::exception &err) {
    AMP::pout << "ERROR: While testing "<<argv[0] << err.what() << std::endl;
    ut.failure("ERROR: While testing");
  } catch( ... ) {
    AMP::pout << "ERROR: While testing "<<argv[0] << "An unknown exception was thrown." << std::endl;
    ut.failure("ERROR: While testing");
  }

  ut.report();

  int num_failed = ut.NumFailGlobal();

  //PROFILE_SAVE(inputFile);

  AMP::AMPManager::shutdown();

  return num_failed;


}   



