#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "ampmesh/Mesh.h"
#include "boost/shared_ptr.hpp"
#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorBuilder.h"
#include "operators/FlowFrapconOperator.h"
#include "operators/FlowFrapconJacobian.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/OperatorBuilder.h"
#include "operators/NeutronicsRhs.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/TrilinosMLSolver.h"
#include "solvers/Flow1DSolver.h"



void PelletCladQuasiStaticThermalFlow(AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);

    // Read the input file
    boost::shared_ptr<AMP::InputDatabase>  input_db ( new AMP::InputDatabase ( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile ( input_file , input_db );
    input_db->printClassData(AMP::plog);

    // Get the Mesh database and create the mesh parameters
    boost::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    boost::shared_ptr<AMP::Mesh::MeshParameters> params(new AMP::Mesh::MeshParameters(database));
    params->setComm(globalComm);

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr  manager = AMP::Mesh::Mesh::buildMesh(params);
    AMP::Mesh::Mesh::shared_ptr  meshAdapter1 = manager->Subset( "pellet" );
    AMP::Mesh::Mesh::shared_ptr  meshAdapter2 = manager->Subset( "clad" );

    //--------------------------------------------------
    // Creating the parameters that will form the right-hand side for the thermal calculation.
    //--------------------------------------------------
    AMP_INSIST(input_db->keyExists("PowerNeutronicsOperator"), "Key ''PowerNeutronicsOperator'' is missing!");
    boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("PowerNeutronicsOperator");
    boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
    neutronicsParams->d_Mesh = meshAdapter1;

    //--------------------------------------------------
    //  Creating a time-steps loop that will be used for the burnup loop.
    //--------------------------------------------------
    int numTimeSteps=3;
    std::vector<double> d_Heff;
    d_Heff.resize(numTimeSteps);
    d_Heff=input_db->getDatabase ( "BoundaryConditions" )->getDoubleArray ( "Heff" );

    //----------------------------------------------------------------------------
    //  Constructing the neutornicsRHS for the Thermal diffusion source (aka specific power).
    //----------------------------------------------------------------------------
    boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator(new AMP::Operator::NeutronicsRhs( neutronicsParams ));

    AMP::LinearAlgebra::Variable::shared_ptr specificPowerGpVar   = neutronicsOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr   specificPowerGpVec   = meshAdapter1->createVector( specificPowerGpVar );

      //----------------------------------------------------------------------------
      //  Create a global temperature variable and an associated mesh specific temperature variable.
      //----------------------------------------------------------------------------
      //AMP::LinearAlgebra::Variable::shared_ptr GlobalTemperatureVar ( new AMP::Mesh::NodalScalarVariable ( "Temperature" ) );
      AMP::LinearAlgebra::Variable::shared_ptr inputThermalVar1 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter1 ) );
      AMP::LinearAlgebra::Variable::shared_ptr inputThermalVar2 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

      AMP::LinearAlgebra::Variable::shared_ptr outputThermalVar1 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter1 ) );
      AMP::LinearAlgebra::Variable::shared_ptr outputThermalVar2 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

      AMP::LinearAlgebra::Variable::shared_ptr thermalMapVar       ( new AMP::Mesh::NodalScalarVariable ( "Temperature" ) );
      AMP::LinearAlgebra::Variable::shared_ptr inputThermalMapVar1 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter1 ) );
      AMP::LinearAlgebra::Variable::shared_ptr inputThermalMapVar2 ( new AMP::Mesh::NodalScalarVariable ( "Temperature", meshAdapter2 ) );

      AMP::LinearAlgebra::Vector::shared_ptr thermalMapVec  = manager->createVector ( thermalMapVar );
      AMP::LinearAlgebra::Vector::shared_ptr thermalMapCladToPelletVec  = thermalMapVec->subsetVectorForVariable ( inputThermalMapVar1  );
      AMP::LinearAlgebra::Vector::shared_ptr thermalMapToCladVec        = thermalMapVec->subsetVectorForVariable ( inputThermalMapVar2  );
      //----------------------------------------------------------------------------
      //  Create a global multivariable with the temperature and displacement on each mesh.
      //----------------------------------------------------------------------------
      boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> globalMultiVar(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
      globalMultiVar->add(inputThermalVar1);
      globalMultiVar->add(inputThermalVar2);

      //----------------------------------------------------------------------------
      //  Create a global multivector with the temperature and displacement on each mesh
      //    for the solution, residual, and right hand side.
      //----------------------------------------------------------------------------
      AMP::LinearAlgebra::Vector::shared_ptr globalSolVec = manager->createVector ( globalMultiVar );
      AMP::LinearAlgebra::Vector::shared_ptr globalRhsVec = manager->createVector ( globalMultiVar );
      AMP::LinearAlgebra::Vector::shared_ptr globalResVec = manager->createVector ( globalMultiVar );

      //-----------------------------------------------
      //   CREATE THE NONLINEAR THERMAL OPERATOR 1 ----
      //-----------------------------------------------
      AMP_INSIST( input_db->keyExists("NonlinearThermalOperator1"), "key missing!" );
      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator1 = boost::dynamic_pointer_cast<
        AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1,
											    "NonlinearThermalOperator1",
											    input_db,
											    thermalTransportModel1));

      AMP::LinearAlgebra::Vector::shared_ptr thermalSolVec1 = globalSolVec->subsetVectorForVariable (  inputThermalVar1 );
      AMP::LinearAlgebra::Vector::shared_ptr thermalRhsVec1 = globalRhsVec->subsetVectorForVariable ( outputThermalVar1 );

      //----------------------------------------------------------------------------
      //  Set the initial guess for the temperature to be 400, or as defined on the input.
      //----------------------------------------------------------------------------
      double intguess = input_db->getDoubleWithDefault("InitialGuess",300);

      thermalSolVec1->setToScalar ( intguess );

      //-------------------------------------
      //   CREATE THE LINEAR THERMAL OPERATOR 1 ----
      //-------------------------------------

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
      boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator1 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1,
																								"LinearThermalOperator1",
																								input_db,
																								thermalTransportModel1));

      //-------------------------------------
      //  CREATE THE NEUTRONICS SOURCE  //
      //-------------------------------------
      AMP::LinearAlgebra::Vector::shared_ptr   nullVec;

      //----------------------------------------------------------
      //  Integrate Nuclear Rhs over Density * Volume //
      //----------------------------------------------------------
      AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );
      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
      boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> specificPowerGpVecToPowerDensityNodalVecOperatator =
        boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1,
															  "VolumeIntegralOperator",
															  input_db,
															  stransportModel));

      //--------------------------------------------
      //   CREATE THE NONLINEAR THERMAL OPERATOR 2 ----
      //--------------------------------------------

      AMP_INSIST( input_db->keyExists("NonlinearThermalOperator2"), "key missing!" );

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
      boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator2 = boost::dynamic_pointer_cast<
        AMP::Operator::NonlinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
											    "NonlinearThermalOperator2",
											    input_db,
											    thermalTransportModel2));

      // initialize the output variable
      // AMP::LinearAlgebra::Variable::shared_ptr outputThermalVariable2 = thermalVolumeOperator2->getOutputVariable();

      AMP::LinearAlgebra::Vector::shared_ptr thermalSolVec2 = globalSolVec->subsetVectorForVariable ( inputThermalVar2  );
      AMP::LinearAlgebra::Vector::shared_ptr thermalRhsVec2 = globalRhsVec->subsetVectorForVariable    ( outputThermalVar2 );

      thermalSolVec2->setToScalar ( intguess );
      //--------------------------------------------


      //--------------------------------------------
      //   CREATE THE LINEAR THERMAL OPERATOR 2 ----
      //--------------------------------------------

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel2;
      boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator2 = boost::dynamic_pointer_cast<
        AMP::Operator::LinearBVPOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
											 "LinearThermalOperator2",
											 input_db,
											 transportModel2));

      //--------------------------------------
      //     CREATE THE FLOW OPERATOR   ------
      //--------------------------------------

      AMP_INSIST(input_db->keyExists("FlowFrapconOperator"), "Key ''FlowFrapconOperator'' is missing!");

      boost::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
      boost::shared_ptr<AMP::InputDatabase> flowDatabase       = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("FlowFrapconOperator"));
      boost::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
																							  "FlowFrapconOperator",
																							  input_db,
																							  flowtransportModel));

      boost::shared_ptr<AMP::Operator::FlowFrapconJacobian> flowJacobian = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
																							  "FlowFrapconJacobian",
																							  input_db,
																							  flowtransportModel));

      //double Cp, De, G, K, Re, Pr, heff, nP;
      double De, K, Re, Pr, heff;

      //Cp  = (flowDatabase)->getDouble("Heat_Capacity");
      De  = (flowDatabase)->getDouble("Channel_Diameter");
      //G   = (flowDatabase)->getDouble("Mass_Flux");
      K   = (flowDatabase)->getDouble("Conductivity");
      Re  = (flowDatabase)->getDouble("Reynolds");
      Pr  = (flowDatabase)->getDouble("Prandtl");
      //nP  = (flowDatabase)->getDouble("numpoints");

      heff = (0.023*K/De)*pow(Re,0.8)*pow(Pr,0.4);


      std::cout << "The flow Heff : "<< heff << std::endl;


      AMP::LinearAlgebra::Vector::shared_ptr flowSol1DVec ;
      //-------------------------------------
      // CREATE MAP OPERATORS
      //-------------------------------------
      boost::shared_ptr<AMP::InputDatabase> mapCladToPellet_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladtoPellet"));
      boost::shared_ptr<AMP::Operator::MapSurface>    mapCladToPellet    = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, meshAdapter1, mapCladToPellet_db ));

      boost::shared_ptr<AMP::InputDatabase> mapPelletToClad_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapPellettoClad"));
      boost::shared_ptr<AMP::Operator::MapSurface>    mapPelletToClad    = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter1, meshAdapter2, mapPelletToClad_db ));

      //---------------------------------------------------------------------------
      //
      //
      //---------------------------------------------------------------------------

      boost::shared_ptr<AMP::InputDatabase> mapcladflow_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladto1DFlow"));
      boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapcladflowParams (new AMP::Operator::MapOperatorParameters( mapcladflow_db ));
      mapcladflowParams->d_MeshAdapter = meshAdapter2;
      boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladTo1DFlow1 (new AMP::Operator::Map3Dto1D( mapcladflowParams ));
      boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladTo1DFlow2 (new AMP::Operator::Map3Dto1D( mapcladflowParams ));

      boost::shared_ptr<AMP::InputDatabase> mapflowclad_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DFlowto3DFlow"));
      boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowcladParams (new AMP::Operator::MapOperatorParameters( mapflowclad_db ));
      mapflowcladParams->d_MapAdapter = meshAdapter2;
      boost::shared_ptr<AMP::Operator::Map1Dto3D> map1DFlowTo3DFlow1 (new AMP::Operator::Map1Dto3D( mapflowcladParams ));
      boost::shared_ptr<AMP::Operator::Map1Dto3D> map1DFlowTo3DFlow2 (new AMP::Operator::Map1Dto3D( mapflowcladParams ));

      mapCladTo1DFlow1->setZLocations( map1DFlowTo3DFlow1->getZLocations());
      mapCladTo1DFlow2->setZLocations( map1DFlowTo3DFlow2->getZLocations());
      size_t flowVecSize = map1DFlowTo3DFlow1->getNumZlocations();

      boost::shared_ptr<AMP::InputDatabase> mapFlowToClad_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapFlowtoClad"));
      boost::shared_ptr<AMP::Operator::MapSurface>    map3DFlowToClad    = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, meshAdapter2, mapFlowToClad_db ));
      //------------------------------------------

      AMP::LinearAlgebra::Vector::shared_ptr cladVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , mapCladTo1DFlow1->getOutputVariable() );
      flowSol1DVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , mapCladTo1DFlow1->getOutputVariable() );

      AMP::LinearAlgebra::Variable::shared_ptr flowVariable = map1DFlowTo3DFlow1->getOutputVariable() ;

      AMP::LinearAlgebra::Vector::shared_ptr flowSolVec = manager->createVector ( flowVariable );
      AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec = manager->createVector ( flowVariable );
      AMP::LinearAlgebra::Vector::shared_ptr flowResVec = manager->createVector ( flowVariable );

      flowOperator->setZLocations( map1DFlowTo3DFlow1->getZLocations());
      flowJacobian->setZLocations( map1DFlowTo3DFlow1->getZLocations());

      flowOperator->setVector(cladVec);
      flowJacobian->setVector(cladVec);

      flowSolVec->setToScalar(300.0);

      //------------------------------------------

      mapCladToPellet->setVector(thermalMapCladToPelletVec );
      mapPelletToClad->setVector(thermalMapToCladVec );
      map3DFlowToClad->setVector(thermalMapToCladVec); 

      mapCladTo1DFlow1->setVector(cladVec); 
      mapCladTo1DFlow2->setVector(cladVec); 
      //------------------------------------------

      boost::shared_ptr<AMP::InputDatabase> tmp1_db (new AMP::InputDatabase("Dummy"));
      boost::shared_ptr<AMP::Operator::OperatorParameters> columnMapsParams (new AMP::Operator::OperatorParameters(tmp1_db));
      boost::shared_ptr<AMP::Operator::ColumnOperator>     columnMapstoCladOperator(new AMP::Operator::ColumnOperator(columnMapsParams));
      columnMapstoCladOperator->append(mapPelletToClad);
      columnMapstoCladOperator->append(map3DFlowToClad);

      //------------------------------------------

      boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp1 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( thermalNonlinearOperator1->getBoundaryOperator()  );

      boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp1->getParameters()) ;

      //------------------------------------------

      boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp2 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( thermalNonlinearOperator2->getBoundaryOperator() ) )->getBoundaryOperator(0) );

      boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters2 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp2->getParameters()) ;

      //------------------------------------------

      boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp3 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( thermalNonlinearOperator2->getBoundaryOperator() ) )->getBoundaryOperator(1) );

      boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters3 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp3->getParameters()) ;

      //------------------------------------------

      robinBoundaryOp1->setVariableFlux( thermalMapCladToPelletVec );
      robinBoundaryOp2->setVariableFlux( thermalMapToCladVec );
      robinBoundaryOp3->setVariableFlux( thermalMapToCladVec );

      robinBoundaryOp1->reset(correctionParameters1);
      robinBoundaryOp2->reset(correctionParameters2);
      robinBoundaryOp3->reset(correctionParameters3);

      //--------------------------------------

      AMP_INSIST(input_db->keyExists("NonlinearSolver"),   "Key ''NonlinearSolver'' is missing!");

      //-------------------------------------
      //Coupling Map to the Nonlinear Operators

      boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
      boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams1(new AMP::Operator::CoupledOperatorParameters(tmp_db));

      coupledNonlinearParams1->d_MapOperator = mapCladToPellet;
      coupledNonlinearParams1->d_BVPOperator = thermalNonlinearOperator1;
      boost::shared_ptr<AMP::Operator::CoupledOperator>           coupledNonlinearOperator1(new AMP::Operator::CoupledOperator(coupledNonlinearParams1));

      //-------------------------------------
      boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams2(new AMP::Operator::CoupledOperatorParameters(tmp_db));
      coupledNonlinearParams2->d_MapOperator = columnMapstoCladOperator;
      coupledNonlinearParams2->d_BVPOperator = thermalNonlinearOperator2;
      boost::shared_ptr<AMP::Operator::CoupledOperator>           coupledNonlinearOperator2(new AMP::Operator::CoupledOperator(coupledNonlinearParams2));

      //-------------------------------------

      boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperatorParameters> coupledNonlinearParams3(new AMP::Operator::CoupledFlowFrapconOperatorParameters(tmp_db));
      coupledNonlinearParams3->d_Map3to1 = mapCladTo1DFlow1 ;
      coupledNonlinearParams3->d_FlowOperator = flowOperator;
      coupledNonlinearParams3->d_Map1to3 = map1DFlowTo3DFlow1;
      coupledNonlinearParams3->d_MeshAdapter = meshAdapter2;
      boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperator>           coupledNonlinearOperator3(new AMP::Operator::CoupledFlowFrapconOperator(coupledNonlinearParams3));

      //-------------------------------------
      //Column of Coupled Operators
      boost::shared_ptr<AMP::Operator::OperatorParameters> columnNonlinearParams (new AMP::Operator::OperatorParameters(tmp_db));
      boost::shared_ptr<AMP::Operator::ColumnOperator>     columnNonlinearOperator(new AMP::Operator::ColumnOperator(columnNonlinearParams));
      columnNonlinearOperator->append(coupledNonlinearOperator1);
      columnNonlinearOperator->append(coupledNonlinearOperator2);
      columnNonlinearOperator->append(coupledNonlinearOperator3);

      //---------------------------------------
      //Column of Coupled Operators
      boost::shared_ptr<AMP::Operator::OperatorParameters> columnLinearParams(new AMP::Operator::OperatorParameters(tmp_db));
      boost::shared_ptr<AMP::Operator::ColumnOperator>     coupledLinearOperator(new AMP::Operator::ColumnOperator(columnLinearParams));
      coupledLinearOperator->append(thermalLinearOperator1);
      coupledLinearOperator->append(thermalLinearOperator2);

      boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperatorParameters> coupledlinearParams3(new AMP::Operator::CoupledFlowFrapconOperatorParameters(tmp_db));
      coupledlinearParams3->d_Map3to1 = mapCladTo1DFlow2 ;
      coupledlinearParams3->d_FlowOperator = flowJacobian;
      coupledlinearParams3->d_Map1to3 = map1DFlowTo3DFlow2;
      coupledlinearParams3->d_MeshAdapter = meshAdapter2;
      boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperator>           coupledlinearOperator3(new AMP::Operator::CoupledFlowFrapconOperator(coupledlinearParams3));

      coupledLinearOperator->append(coupledlinearOperator3);

      //----------------------------------------------------------------------------------------------------------------------------------------------//

      AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
      globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalSolVec );
      globalSolMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowSolVec );

      AMP::LinearAlgebra::Vector::shared_ptr globalSolMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalSolMultiVector, globalComm );
      //---------------------------------------------------------------------------------------------------------------------//
      AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
      globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalRhsVec );
      globalRhsMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowRhsVec );

      AMP::LinearAlgebra::Vector::shared_ptr globalRhsMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalRhsMultiVector, globalComm );
      //---------------------------------------------------------------------------------------------------------------------//
      //---------------------------------------------------------------------------------------------------------------------//
      AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVector = AMP::LinearAlgebra::MultiVector::create( "multivector" , globalComm ) ;
      globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( globalResVec );
      globalResMultiVector->castTo<AMP::LinearAlgebra::MultiVector>().addVector ( flowResVec );

//      AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalResMultiVector, globalComm );
      //---------------------------------------------------------------------------------------------------------------------//
      //---------------------------------------------------------------------------------------------------------------------//

      neutronicsOperator->setTimeStep(0);
      neutronicsOperator->apply(nullVec, nullVec, specificPowerGpVec, 1., 0.);
      thermalRhsVec1->zero();
      specificPowerGpVecToPowerDensityNodalVecOperatator->apply(nullVec, specificPowerGpVec, thermalRhsVec1, 1., 0.);

      //We need to reset the linear operator before the solve since TrilinosML does
      //the factorization of the matrix during construction and so the matrix must
      //be correct before constructing the TrilinosML object.
      //The thermal operator does not expect an apply to be called before calling
      //getJacobianParams and so it need not be called. So, any of the following
      //apply calls will work:
      coupledLinearOperator->reset(columnNonlinearOperator->getJacobianParameters(globalSolMultiVector));
      columnNonlinearOperator->apply(nullVec, globalSolMultiVector, globalResMultiVector, 1.0, 0.0);

      //------------------------------------------------------------------
      boost::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase("NonlinearSolver");
      boost::shared_ptr<AMP::Database>    linearSolver_db = nonlinearSolver_db->getDatabase("LinearSolver");

      //----------------------------------------------------------------//
      // initialize the nonlinear solver
      boost::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(new AMP::Solver::PetscSNESSolverParameters(nonlinearSolver_db));

      // change the next line to get the correct communicator out
      nonlinearSolverParams->d_comm          = globalComm;
      nonlinearSolverParams->d_pOperator     = columnNonlinearOperator;
      nonlinearSolverParams->d_pInitialGuess = globalSolMultiVector ;
      boost::shared_ptr<AMP::Solver::PetscSNESSolver>  nonlinearSolver(new AMP::Solver::PetscSNESSolver(nonlinearSolverParams));

      //-------------------------------------------------------------------------//
      // initialize the column preconditioner which is a diagonal block preconditioner
      boost::shared_ptr<AMP::Database>                 columnPreconditioner_db = linearSolver_db->getDatabase("Preconditioner");
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> columnPreconditionerParams(new AMP::Solver::SolverStrategyParameters(columnPreconditioner_db));
      columnPreconditionerParams->d_pOperator = coupledLinearOperator;
      boost::shared_ptr<AMP::Solver::ColumnSolver>             columnPreconditioner(new AMP::Solver::ColumnSolver(columnPreconditionerParams));

      //---------------
      boost::shared_ptr<AMP::Database>                 thermalPreconditioner_db1 = columnPreconditioner_db->getDatabase("pelletThermalPreconditioner");
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams1(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db1));
      thermalPreconditionerParams1->d_pOperator = thermalLinearOperator1;
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         thermalPreconditioner1(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams1));

      boost::shared_ptr<AMP::Database>                 thermalPreconditioner_db2 = columnPreconditioner_db->getDatabase("cladThermalPreconditioner");
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams2(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db2));
      thermalPreconditionerParams2->d_pOperator = thermalLinearOperator2;
      boost::shared_ptr<AMP::Solver::TrilinosMLSolver>         thermalPreconditioner2(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams2));
      //---------------
      boost::shared_ptr<AMP::Database> JacobianSolver_db = input_db->getDatabase("Flow1DSolver"); 
      boost::shared_ptr<AMP::Solver::SolverStrategyParameters> flowSolverParams (new AMP::Solver::SolverStrategyParameters(JacobianSolver_db));
      flowSolverParams->d_pOperator = flowJacobian;
      boost::shared_ptr<AMP::Solver::Flow1DSolver>         flowJacobianSolver(new AMP::Solver::Flow1DSolver(flowSolverParams));

      boost::shared_ptr<AMP::InputDatabase> CoupledJacobianSolver_db (new AMP::InputDatabase("Dummy"));
      CoupledJacobianSolver_db->putInteger("max_iterations",1); 
      CoupledJacobianSolver_db->putDouble("max_error",1.0e-6); 
      boost::shared_ptr<AMP::Solver::CoupledFlow1DSolverParameters> coupledFlowSolverParams (new AMP::Solver::CoupledFlow1DSolverParameters(CoupledJacobianSolver_db ));
      coupledFlowSolverParams->d_flow1DSolver = flowJacobianSolver;
      coupledFlowSolverParams->d_pOperator = boost::dynamic_pointer_cast< AMP::Operator::Operator>(coupledlinearOperator3);
      boost::shared_ptr<AMP::Solver::CoupledFlow1DSolver>   CoupledFlowJacobianSolver(new AMP::Solver::CoupledFlow1DSolver(coupledFlowSolverParams));
      //---------------

      columnPreconditioner->append(thermalPreconditioner1);
      columnPreconditioner->append(thermalPreconditioner2);
      columnPreconditioner->append(CoupledFlowJacobianSolver);

      //--------------------------------------------------------------------//
      // register the preconditioner with the Jacobian free Krylov solver
      boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
      linearSolver->setPreconditioner(columnPreconditioner);

      //-------------------------------------
      nonlinearSolver->setZeroInitialGuess(false);

      AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  b = meshAdapter2->beginOwnedBoundary ( 4 );
      AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator  e = meshAdapter2->endOwnedBoundary ( 4 );
      //AMP_ASSERT ( b != e );
      AMP::LinearAlgebra::Vector::shared_ptr surfaceVector = flowSolVec->select ( AMP::Mesh::VS_ByMeshIteratorTmpl< AMP::Mesh::MeshManager::Adapter , AMP::Mesh::MeshManager::Adapter::OwnedBoundaryNodeIterator >
( meshAdapter2 , b , e , 1 ) , "Flow" );


      manager->registerVectorAsData ( thermalSolVec1 , "PelletTemperature");
      manager->registerVectorAsData ( thermalSolVec2 , "CladTemperature");
      manager->registerVectorAsData ( surfaceVector ,     "FLowTemperature");

      for ( int tstep = 0; tstep < 1; tstep++ ) {
        neutronicsOperator->setTimeStep(tstep);
        neutronicsOperator->apply(nullVec, nullVec, specificPowerGpVec, 1., 0.);

        robinBoundaryOp1->reset(correctionParameters1);
        robinBoundaryOp2->reset(correctionParameters2);
        robinBoundaryOp3->reset(correctionParameters3);

        thermalRhsVec1->zero();
        // specificPowerGpVec is in Watts/kilogram
        specificPowerGpVecToPowerDensityNodalVecOperatator->apply(nullVec, specificPowerGpVec, thermalRhsVec1, 1., 0.);

        thermalNonlinearOperator1->modifyRHSvector(thermalRhsVec1);
        thermalNonlinearOperator1->modifyInitialSolutionVector(thermalSolVec1);
        thermalNonlinearOperator2->modifyRHSvector(thermalRhsVec2);
        thermalNonlinearOperator2->modifyInitialSolutionVector(thermalSolVec2);

        columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
        AMP::pout<<"Initial Residual Norm for Step " << tstep << " is: "<<globalResMultiVector->L2Norm()<<std::endl;

        nonlinearSolver->solve(globalRhsMultiVectorView, globalSolMultiVectorView);

        columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
        AMP::pout<<"Final   Residual Norm for Step " << tstep << " is: "<<globalResMultiVector->L2Norm()<<std::endl;
#ifdef USE_SILO
        manager->writeFile<AMP::Mesh::SiloIO> ( silo_name , tstep );
#endif

        std::cout<< "The Fuel/clad Max:Min values - "<<thermalSolVec1->max() << " " <<thermalSolVec2->min() <<std::endl;
        std::cout << "Surface Vector Max:Min values -  " << surfaceVector->max() << " " << surfaceVector->min() <<std::endl;

        std::cout<<"Intermediate Flow Solution " <<std::endl;
        mapCladTo1DFlow1->setVector(flowSol1DVec); 
	mapCladTo1DFlow1->apply(nullVec, thermalMapToCladVec , nullVec, 1, 0);
        std::vector<double> expectedSolution(flowVecSize, 0);
        expectedSolution = (input_db->getDatabase("regression"))->getDoubleArray("expectedSolution");
        for(unsigned int i=0; i<flowVecSize ; i++) {
          if( !AMP::Utilities::approx_equal(expectedSolution[i], flowSol1DVec->getValueByLocalID(i), 1e-6)) {

            if(globalComm.getRank() == 0){
              printf("solution: %.7e expected: %.7e \n", 
                flowSol1DVec->getValueByLocalID(i),expectedSolution[i] );
            
            }
            ut->failure("solution is different for "+silo_name);
          }
        }
        std::cout<<std::endl;

        nonlinearSolver->setZeroInitialGuess(false);
      }

      input_db.reset();

      if( ut->NumFailLocal() ==0 ) ut->passes(exeName);

}

int main(int argc, char *argv[])
{
  AMP::AMPManager::startup(argc, argv);
  AMP::UnitTest ut;

  try {
    PelletCladQuasiStaticThermalFlow(&ut, "testFlowSolution-2");
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


