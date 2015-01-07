#include <utils/InputDatabase.h>
#include <utils/InputManager.h>
#include <utils/AMP_MPI.h>
#include <utils/AMPManager.h>
#include <utils/UnitTest.h>
#include <utils/Writer.h>
#include <utils/Utilities.h>
#include <utils/PIO.h>
#include <ampmesh/Mesh.h>
#include <discretization/simpleDOF_Manager.h>
#include <operators/map/dtk/DTKMapOperator.h>
#include <operators/ColumnOperator.h>
#include <vectors/Variable.h>
#include <vectors/VectorBuilder.h>
#include <vectors/MultiVector.h>
#include <vectors/VectorSelector.h>

#define __INIT_FN1__(x, y, z) ( x+y+z )
#define __INIT_FN2__(x, y, z) ( 2*x+y+z )
#define __INIT_FN3__(x, y, z) ( 4*x+y+z )

int runTest(std::string exeName, AMP::UnitTest *ut)
{
    std::string const inputFile = "input_" + exeName;
    std::string const logFile   = "output_" + exeName;

    AMP::PIO::logAllNodes(logFile);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Parse input file
    AMP::shared_ptr<AMP::InputDatabase> inputDatabase(new AMP::InputDatabase("input_db"));
    AMP::InputManager::getManager()->parseInputFile(inputFile, inputDatabase);
    inputDatabase->printClassData(AMP::pout);

    // Read the mesh
    AMP::pout<<"--------------------\n";
    AMP::pout<<"    LOADING MESH    \n";
    AMP::pout<<"--------------------\n";
    AMP::Database::shared_ptr meshDatabase = inputDatabase->getDatabase("Mesh");
    AMP::Mesh::MeshParameters::shared_ptr meshParams(new AMP::Mesh::MeshParameters(meshDatabase));
    meshParams->setComm(globalComm);
    AMP::Mesh::Mesh::shared_ptr mesh = AMP::Mesh::Mesh::buildMesh(meshParams);

    // Subset the mesh
    AMP::Mesh::Mesh::shared_ptr cellSandwichMesh = mesh->Subset("CellSandwich_2_1");
    AMP::Mesh::Mesh::shared_ptr anodeCCMesh      = mesh->Subset("AnodeCC_1_1"     );
    AMP::Mesh::Mesh::shared_ptr cathodeCCMesh    = mesh->Subset("CathodeCC_3_1"   );

    // Distribute degrees of freedom
    int const ghostWidth = 0;
    bool const split = true;
    AMP::Discretization::DOFManager::shared_ptr phiDofMap  = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, 1, split);
    AMP::Discretization::DOFManager::shared_ptr eectDofMap = 
        AMP::Discretization::simpleDOFManager::create(mesh, AMP::Mesh::Vertex, ghostWidth, 5, split);

  // Construct Vectors
  AMP::pout<<"------------------------------------\n";
  AMP::pout<<"     BUILD VECTORS AND FILL THEM    \n";
  AMP::pout<<"------------------------------------\n";
  AMP::LinearAlgebra::Variable::shared_ptr potentialVariable(new AMP::LinearAlgebra::Variable("Potential"));
  AMP::LinearAlgebra::Vector::shared_ptr   potentialMapVec = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, split);
  AMP::LinearAlgebra::Vector::shared_ptr   potentialSolVec = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, split);
  AMP::LinearAlgebra::Vector::shared_ptr   potentialResVec = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, split);
  AMP::LinearAlgebra::Vector::shared_ptr   potentialRhsVec = AMP::LinearAlgebra::createVector( phiDofMap , potentialVariable, split);

  AMP::LinearAlgebra::Variable::shared_ptr batteryVariables (new AMP::LinearAlgebra::Variable("Battery"));
  AMP::LinearAlgebra::Vector::shared_ptr   BatterySolVec   = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , split);
  AMP::LinearAlgebra::Vector::shared_ptr   BatteryResVec   = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , split);
  AMP::LinearAlgebra::Vector::shared_ptr   BatteryMapVec   = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , split);
  AMP::LinearAlgebra::Vector::shared_ptr   BatteryRhsVec   = AMP::LinearAlgebra::createVector( eectDofMap    , batteryVariables , split);

  AMP::LinearAlgebra::Vector::shared_ptr   ElectrodeSolVec = BatterySolVec->select( AMP::LinearAlgebra::VS_Stride( 3, 5) , "V4" );
  AMP::LinearAlgebra::Vector::shared_ptr   ElectrodeMapVec = BatteryMapVec->select( AMP::LinearAlgebra::VS_Stride( 3, 5) , "V4" );
//---------------------------------------------------

  AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
//  siloWriter->registerMesh( mesh );
  siloWriter->setDecomposition(1);
  siloWriter->registerVector(potentialMapVec, mesh, AMP::Mesh::Vertex, "potentialMapVec");
  siloWriter->registerVector(potentialSolVec, mesh, AMP::Mesh::Vertex, "potentialSolVec");
  siloWriter->registerVector(ElectrodeMapVec, mesh, AMP::Mesh::Vertex, "batteryMapVec"  );
  siloWriter->registerVector(ElectrodeSolVec, mesh, AMP::Mesh::Vertex, "batterySolVec"  );
  siloWriter->writeFile(logFile , 0);

//---------------------------------------------------

  AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiSolVec = AMP::LinearAlgebra::MultiVector::create("MultiSolVec", globalComm);
  multiSolVec->addVector(BatterySolVec);
  multiSolVec->addVector(potentialSolVec);

  AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiResVec = AMP::LinearAlgebra::MultiVector::create("MultiResVec", globalComm);
  multiResVec->addVector(BatteryResVec);
  multiResVec->addVector(potentialResVec);

  AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiRhsVec = AMP::LinearAlgebra::MultiVector::create("MultiRhsVec", globalComm);
  multiRhsVec->addVector(BatteryRhsVec);
  multiRhsVec->addVector(potentialRhsVec);

   // Make new vectors
  AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> multiMapVec = AMP::LinearAlgebra::MultiVector::create("MultiMapVec", globalComm);
  multiMapVec->addVector(BatteryMapVec);
  multiMapVec->addVector(potentialMapVec);



  // Filling the vectors
  //     potential on the electrodes
  AMP::Mesh::MeshIterator node  = anodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 5, 0 );
  AMP::Mesh::MeshIterator end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN1__(pt[0], pt[1], pt[2]);
        AMP_ASSERT(bndGlobalIds.size() == 1);

        potentialSolVec->setValueByGlobalID(bndGlobalIds[0], val);
  }//end for node

  node  = cathodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 3, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN3__(pt[0], pt[1], pt[2]);
        AMP_ASSERT(bndGlobalIds.size() == 1);

        potentialSolVec->setValueByGlobalID(bndGlobalIds[0], val);
  }//end for node

  //     4th component on the sandwich
  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN2__(pt[0], pt[1], pt[2]);
        AMP_ASSERT(bndGlobalIds.size() == 5);

        BatterySolVec->setValueByGlobalID(bndGlobalIds[3], val);
  }//end for node

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
  end_node = node.end();
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN2__(pt[0], pt[1], pt[2]);
        AMP_ASSERT(bndGlobalIds.size() == 5);

        BatterySolVec->setValueByGlobalID(bndGlobalIds[3], val);
  }//end for node
  multiSolVec->makeConsistent(AMP::LinearAlgebra::Vector::CONSISTENT_SET);

  
  siloWriter->writeFile( logFile , 1);

    // create dtk map operator.
    AMP::pout<<"----------------------------\n";
    AMP::pout<<"     CREATE MAP OPERATOR    \n";
    AMP::pout<<"----------------------------\n";

    ///////////////////////////////////////////////////////////////////
    // TODO: MAKE AN OPERATOR THAT TAKE THIS AS PARAMETERS
    AMP::Mesh::Mesh::shared_ptr mesh1 = cellSandwichMesh;
    AMP::Mesh::Mesh::shared_ptr mesh2 = cathodeCCMesh   ;
    AMP::LinearAlgebra::Vector::shared_ptr sourceVector = multiSolVec;
    AMP::LinearAlgebra::Vector::shared_ptr targetVector = multiMapVec;
    int const boundaryID1 = 2;
    int const boundaryID2 = 3;
    std::string const variable1 = batteryVariables ->getName();
    std::size_t const strideOffset1 = 3;
    std::size_t const strideLength1 = 5;
    std::string const variable2 = potentialVariable->getName();
    std::size_t const strideOffset2 = 0;
    std::size_t const strideLength2 = 1;
    ///////////////////////////////////////////////////////////////////


    AMP::shared_ptr<AMP::Database> nullDatabase;
    // map 1 -> 2
    AMP::LinearAlgebra::Vector::shared_ptr sourceVectorMap12 = sourceVector
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable1)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset1, strideLength1), "var");
    AMP::LinearAlgebra::Vector::shared_ptr targetVectorMap12 = targetVector
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable1)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset1, strideLength1), "var");
    AMP::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map12Params(new AMP::Operator::DTKMapOperatorParameters(nullDatabase));
    map12Params->d_domain_mesh = mesh1->Subset(mesh1->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID1));
    map12Params->d_range_mesh  = mesh2->Subset(mesh2->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID2));
    map12Params->d_domain_dofs = sourceVectorMap12->getDOFManager();
    map12Params->d_range_dofs  = targetVectorMap12->getDOFManager();
    AMP::shared_ptr<AMP::Operator::Operator> map12(new AMP::Operator::DTKMapOperator(map12Params));

    // map 2 -> 1
    AMP::LinearAlgebra::Vector::shared_ptr sourceVectorMap21 = sourceVector
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable2)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset2, strideLength2), "var");
    AMP::LinearAlgebra::Vector::shared_ptr targetVectorMap21 = targetVector
            ->select(AMP::LinearAlgebra::VS_ByVariableName(variable2)           , "var")
            ->select(AMP::LinearAlgebra::VS_Stride(strideOffset2, strideLength2), "var");
    AMP::shared_ptr<AMP::Operator::DTKMapOperatorParameters> map21Params(new AMP::Operator::DTKMapOperatorParameters(nullDatabase));
    map21Params->d_domain_mesh = mesh2->Subset(mesh2->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID2));
    map21Params->d_range_mesh  = mesh1->Subset(mesh1->getBoundaryIDIterator(AMP::Mesh::Volume, boundaryID1));
    map21Params->d_domain_dofs = sourceVectorMap21->getDOFManager();
    map21Params->d_range_dofs  = targetVectorMap21->getDOFManager();
    AMP::shared_ptr<AMP::Operator::Operator> map21(new AMP::Operator::DTKMapOperator(map21Params));

    // apply the map.
    AMP::pout<<"----------------------\n";
    AMP::pout<<"     APPLY THE MAP    \n";
    AMP::pout<<"----------------------\n";

    ////////////////////////////////////////////////////////
    // TODO: DOES TO THE APPLY
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    map12->apply(nullVec, sourceVectorMap12, targetVectorMap12);
    map21->apply(nullVec, sourceVectorMap21, targetVectorMap21);
    ////////////////////////////////////////////////////////

    siloWriter->writeFile( logFile , 2);


    // check the answer
    AMP::pout<<"----------------------\n";
    AMP::pout<<"     COMPUTE ERROR    \n";
    AMP::pout<<"----------------------\n";
    double const absoluteTolerance = 1.0e-14;
    double const relativeTolerance = 1.0e-14;
    AMP::LinearAlgebra::Vector::const_shared_ptr errorVec;
    double tolerance;
    double errorNorm;
    std::string whatAmIChecking;

    // FIRST
    node  = cathodeCCMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 3, 0 );
    end_node = node.end();
    errorVec = ElectrodeMapVec->constSelect(AMP::LinearAlgebra::VS_MeshIterator(node.begin(), globalComm), "error");
    tolerance = absoluteTolerance + relativeTolerance * errorVec->L2Norm();
    whatAmIChecking = "interface between cellSandwich and cathodeCC map electrochemical";
    
    
  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        eectDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN2__(pt[0], pt[1], pt[2]);
        BatteryMapVec->addValueByGlobalID(bndGlobalIds[3], -val);
  }//end for node
    errorNorm = errorVec->L2Norm();
    AMP::pout<<"anodeCC error = "<<errorNorm<<"\n"; 
    if (errorNorm > tolerance)
        ut->failure(whatAmIChecking);

    // SECOND
    node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
    end_node = node.end();
    errorVec = potentialMapVec->constSelect(AMP::LinearAlgebra::VS_MeshIterator(node.begin(), globalComm), "error");
    tolerance = absoluteTolerance + relativeTolerance * errorVec->L2Norm();
    whatAmIChecking = "interface between cellSandwich and cathodeCC map electrical";

  for( ; node != end_node ; ++node)
  {
        std::vector<size_t> bndGlobalIds;
        phiDofMap->getDOFs( node->globalID() , bndGlobalIds );
        
        std::vector<double> pt = node->coord();
        double val = __INIT_FN3__(pt[0], pt[1], pt[2]);
        potentialMapVec->addValueByGlobalID(bndGlobalIds[0], -val);
  }//end for node
    errorNorm = errorVec->L2Norm();
    AMP::pout<<"sandwich anodeCC error = "<<errorNorm<<"\n"; 
    if (errorNorm > tolerance)
        ut->failure(whatAmIChecking);

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 1, 0 );
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

  node  = cellSandwichMesh->getBoundaryIDIterator( AMP::Mesh::Vertex, 2, 0 );
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

  siloWriter->writeFile( logFile , 3);

  return 1;
}

int main(int argc, char *argv[])
{
    AMP::AMPManagerProperties startup_properties;
    startup_properties.use_MPI_Abort = false;
    AMP::AMPManager::startup(argc,argv,startup_properties);
    AMP::UnitTest ut;

    runTest("testMultiDofDTKMap-1", &ut);

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
