
#include <string>
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "boost/shared_ptr.hpp"

#include "utils/InputDatabase.h"
#include "utils/Utilities.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/Database.h"

#include "ampmesh/Mesh.h"
#include "utils/Writer.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/VectorBuilder.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorSelector.h"
#include "operators/ColumnOperator.h"
#include "operators/CoupledOperator.h"
#include "operators/subchannel/CoupledFlowFrapconOperator.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/ElementOperationFactory.h"
#include "operators/OperatorBuilder.h"
#include "operators/NeutronicsRhs.h"
#include "operators/LinearBVPOperator.h"
#include "operators/NonlinearBVPOperator.h"
#include "operators/VolumeIntegralOperator.h"
#include "operators/boundary/RobinVectorCorrection.h"
#include "operators/boundary/ColumnBoundaryOperator.h"
#include "operators/map/Map3Dto1D.h"
#include "operators/map/Map1Dto3D.h"
#include "operators/map/AsyncMapColumnOperator.h"
#include "operators/map/ScalarZAxisMap.h"
#include "operators/map/MapSurface.h"
#include "solvers/ColumnSolver.h"
#include "solvers/PetscKrylovSolverParameters.h"
#include "solvers/PetscKrylovSolver.h"
#include "solvers/PetscSNESSolverParameters.h"
#include "solvers/PetscSNESSolver.h"
#include "solvers/trilinos/TrilinosMLSolver.h"
#include "solvers/Flow1DSolver.h"
#include "solvers/CoupledFlow1DSolver.h"


void PelletCladQuasiStaticThermalFlow(AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file = "output_" + exeName;
    std::string silo_name = exeName;
    AMP::PIO::logAllNodes(log_file);
    AMP::AMP_MPI globalComm(AMP_COMM_WORLD);
    globalComm.barrier();
    double t0 = AMP::AMP_MPI::time();

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

    // Create a surface mesh on the clad
    AMP::Mesh::Mesh::shared_ptr surfaceMesh;
    if ( meshAdapter2.get() != NULL ) {
        surfaceMesh = meshAdapter2->Subset( meshAdapter2->getBoundaryIDIterator( AMP::Mesh::Face, 4, 0 ) );
        surfaceMesh->setName( "clad_surface" );
    }
    globalComm.barrier();
    double t1 = AMP::AMP_MPI::time();
    std::cout << "Time to load meshes: " << t1-t0 << std::endl;

    // Create the DOF managers
    AMP::Discretization::DOFManager::shared_ptr nodalScalarDOF = 
        AMP::Discretization::simpleDOFManager::create(manager,AMP::Mesh::Vertex,1,1,true);
    //AMP::Discretization::DOFManager::shared_ptr flowNodalScalarDOF = 
    //    AMP::Discretization::simpleDOFManager::create(surfaceMesh,AMP::Mesh::Vertex,1,1,true);
    AMP::Discretization::DOFManager::shared_ptr flowNodalScalarDOF;
    if ( meshAdapter2.get()!=NULL )
        flowNodalScalarDOF = AMP::Discretization::simpleDOFManager::create(meshAdapter2,AMP::Mesh::Vertex,1,1,true);
    int DOFsPerElement = 8;
    AMP::Discretization::DOFManager::shared_ptr gaussPointDOF1;
    if ( meshAdapter1.get()!=NULL )
        gaussPointDOF1 = AMP::Discretization::simpleDOFManager::create(meshAdapter1,AMP::Mesh::Volume,1,DOFsPerElement,true);

    //--------------------------------------------------
    // Creating the parameters that will form the right-hand side for the thermal calculation.
    //--------------------------------------------------
    AMP_INSIST(input_db->keyExists("PowerNeutronicsOperator"), "Key ''PowerNeutronicsOperator'' is missing!");
    boost::shared_ptr<AMP::Database>  neutronicsOp_db = input_db->getDatabase("PowerNeutronicsOperator");
    boost::shared_ptr<AMP::Operator::NeutronicsRhsParameters> neutronicsParams(new AMP::Operator::NeutronicsRhsParameters( neutronicsOp_db ));
    neutronicsParams->d_Mesh = meshAdapter1;

    //----------------------------------------------------------------------------
    //  Constructing the neutornicsRHS for the Thermal diffusion source (aka specific power).
    //----------------------------------------------------------------------------
    boost::shared_ptr<AMP::Operator::NeutronicsRhs> neutronicsOperator;
    AMP::LinearAlgebra::Vector::shared_ptr   specificPowerGpVec;
    if ( meshAdapter1.get()!=NULL ) {
        neutronicsOperator = boost::shared_ptr<AMP::Operator::NeutronicsRhs>(new AMP::Operator::NeutronicsRhs( neutronicsParams ));
        AMP::LinearAlgebra::Variable::shared_ptr specificPowerGpVar = neutronicsOperator->getOutputVariable();
        specificPowerGpVec = AMP::LinearAlgebra::createVector( gaussPointDOF1, specificPowerGpVar, true );
    }

    //----------------------------------------------------------------------------
    //  Create a global temperature variable and an associated mesh specific temperature variable.
    //----------------------------------------------------------------------------
    //AMP::LinearAlgebra::Variable::shared_ptr GlobalTemperatureVar ( new AMP::Mesh::NodalScalarVariable ( "Temperature" ) );
    AMP::LinearAlgebra::Variable::shared_ptr thermalVar ( new AMP::LinearAlgebra::Variable( "Temperature" ) );

    AMP::LinearAlgebra::Vector::shared_ptr  thermalMapVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, thermalVar, true );

    AMP::LinearAlgebra::Vector::shared_ptr thermalMapCladToPelletVec;
    AMP::LinearAlgebra::Vector::shared_ptr thermalMapToCladVec;
    if ( meshAdapter1.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector1( meshAdapter1 );
        thermalMapCladToPelletVec = thermalMapVec->select( meshSelector1, "Temperature" );
    }
    if ( meshAdapter2.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector2( meshAdapter2 );
        thermalMapToCladVec = thermalMapVec->select( meshSelector2, "Temperature" );
    }

    //----------------------------------------------------------------------------
    //  Create a global multivariable with the temperature and displacement on each mesh.
    //----------------------------------------------------------------------------
    boost::shared_ptr<AMP::LinearAlgebra::MultiVariable> globalMultiVar(new AMP::LinearAlgebra::MultiVariable("inputVariable"));
    globalMultiVar->add(thermalVar);

    //----------------------------------------------------------------------------
    //  Create a global multivector with the temperature and displacement on each mesh
    //    for the solution, residual, and right hand side.
    //----------------------------------------------------------------------------
    AMP::LinearAlgebra::Vector::shared_ptr globalSolVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, globalMultiVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr globalRhsVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, globalMultiVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr globalResVec = AMP::LinearAlgebra::createVector( nodalScalarDOF, globalMultiVar, true );
    AMP::LinearAlgebra::Vector::shared_ptr thermalSolVec1, thermalRhsVec1, thermalSolVec2, thermalRhsVec2;
    double intguess = input_db->getDoubleWithDefault("InitialGuess",300);
    if ( meshAdapter1.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector1( meshAdapter1 );
        thermalSolVec1 = globalSolVec->select( meshSelector1, "Temperature" );
        //  Set the initial guess for the temperature to be 400, or as defined on the input.
        thermalRhsVec1 = globalRhsVec->select( meshSelector1, "Temperature" );
        thermalSolVec1->setToScalar ( intguess );
    }
    if ( meshAdapter2.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector2( meshAdapter2 );
        thermalSolVec2 = globalSolVec->select( meshSelector2, "Temperature" );
        thermalRhsVec2 = globalRhsVec->select( meshSelector2, "Temperature" );
        thermalSolVec2->setToScalar ( intguess );
    }

    //-------------------------------------
    //  CREATE THE NEUTRONICS SOURCE  //
    //-------------------------------------
    AMP::LinearAlgebra::Vector::shared_ptr   nullVec;

    //-----------------------------------------------
    //   CREATE THE NONLINEAR AND LINEAR THERMAL OPERATOR 1 ----
    //-----------------------------------------------
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator1;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator1;
    if ( meshAdapter1.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector1( meshAdapter1 );
        AMP_INSIST( input_db->keyExists("NonlinearThermalOperator1"), "key missing!" );
        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel1;
        thermalNonlinearOperator1 = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( 
            AMP::Operator::OperatorBuilder::createOperator( meshAdapter1, "NonlinearThermalOperator1", input_db, thermalTransportModel1));

        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel1;
        thermalLinearOperator1 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( meshAdapter1, "LinearThermalOperator1", input_db, thermalTransportModel1) );
    }

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Density * Volume //
    //----------------------------------------------------------
    boost::shared_ptr<AMP::Operator::VolumeIntegralOperator> specificPowerGpVecToPowerDensityNodalVecOperatator;
    if ( meshAdapter1.get()!=NULL ) {
        AMP_INSIST( input_db->keyExists("VolumeIntegralOperator"), "key missing!" );
        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
        specificPowerGpVecToPowerDensityNodalVecOperatator = boost::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator( meshAdapter1, "VolumeIntegralOperator", input_db, stransportModel) );
    }

    //--------------------------------------------
    //   CREATE THE NONLINEAR AND LINEAR THERMAL OPERATOR 2 ----
    //--------------------------------------------
    boost::shared_ptr<AMP::Operator::NonlinearBVPOperator> thermalNonlinearOperator2;
    boost::shared_ptr<AMP::Operator::LinearBVPOperator> thermalLinearOperator2;
    if ( meshAdapter2.get()!=NULL ) {
        AMP::LinearAlgebra::VS_Mesh meshSelector2( meshAdapter2 );
        AMP_INSIST( input_db->keyExists("NonlinearThermalOperator2"), "key missing!" );

        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> thermalTransportModel2;
        thermalNonlinearOperator2 = boost::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
											    "NonlinearThermalOperator2",
											    input_db,
											    thermalTransportModel2));

        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel2;
        thermalLinearOperator2 = boost::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( meshAdapter2, "LinearThermalOperator2", input_db, transportModel2 ) );
    }

    //--------------------------------------
    //     CREATE THE FLOW OPERATOR   ------
    //--------------------------------------
    boost::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator;
    boost::shared_ptr<AMP::Operator::FlowFrapconJacobian> flowJacobian;
    if ( meshAdapter2.get()!=NULL ) {
        AMP_INSIST(input_db->keyExists("FlowFrapconOperator"), "Key ''FlowFrapconOperator'' is missing!");

        boost::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
        flowOperator = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
																							  "FlowFrapconOperator",
																							  input_db,
																							  flowtransportModel));

        flowJacobian = boost::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>(AMP::Operator::OperatorBuilder::createOperator(meshAdapter2,
																							  "FlowFrapconJacobian",
																							  input_db,
																							  flowtransportModel));
    }

    boost::shared_ptr<AMP::InputDatabase> flowDatabase = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("FlowFrapconOperator"));

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



    //-------------------------------------
    // CREATE THE PELLET-CLAD MAP OPERATORS
    //-------------------------------------
    boost::shared_ptr<AMP::Operator::AsyncMapColumnOperator>  mapPelletAndClad;
    mapPelletAndClad = AMP::Operator::AsyncMapColumnOperator::build<AMP::Operator::ScalarZAxisMap> ( manager, input_db->getDatabase("PelletAndCladMaps") );
    mapPelletAndClad->setVector( thermalMapVec );


    //---------------------------------------------------------------------------
    // CREATE THE FLOW MAP OPERATORS
    //---------------------------------------------------------------------------
    boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladTo1DFlow1;
    boost::shared_ptr<AMP::Operator::Map3Dto1D> mapCladTo1DFlow2;
    boost::shared_ptr<AMP::Operator::Map1Dto3D> map1DFlowTo3DFlow1;
    boost::shared_ptr<AMP::Operator::Map1Dto3D> map1DFlowTo3DFlow2;
    boost::shared_ptr<AMP::Operator::MapSurface>  map3DFlowToClad;
    AMP::LinearAlgebra::Vector::shared_ptr cladVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowSol1DVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowSolVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowRhsVec;
    AMP::LinearAlgebra::Vector::shared_ptr flowResVec;
    if ( meshAdapter2.get()!=NULL ) {
        boost::shared_ptr<AMP::InputDatabase> mapcladflow_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapCladto1DFlow"));
        boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapcladflowParams (new AMP::Operator::MapOperatorParameters( mapcladflow_db ));
        mapcladflowParams->d_Mesh = meshAdapter2;
        mapcladflowParams->d_MapMesh = meshAdapter2;
        mapcladflowParams->d_MapComm = meshAdapter2->getComm();  
        mapCladTo1DFlow1 = boost::shared_ptr<AMP::Operator::Map3Dto1D>(new AMP::Operator::Map3Dto1D( mapcladflowParams ));
        mapCladTo1DFlow2 = boost::shared_ptr<AMP::Operator::Map3Dto1D>(new AMP::Operator::Map3Dto1D( mapcladflowParams ));

        boost::shared_ptr<AMP::InputDatabase> mapflowclad_db  = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("Map1DFlowto3DFlow"));
        boost::shared_ptr<AMP::Operator::MapOperatorParameters> mapflowcladParams (new AMP::Operator::MapOperatorParameters( mapflowclad_db ));
        mapflowcladParams->d_Mesh = meshAdapter2;
        mapflowcladParams->d_MapMesh = meshAdapter2;
        mapflowcladParams->d_MapComm = meshAdapter2->getComm(); 
        map1DFlowTo3DFlow1 = boost::shared_ptr<AMP::Operator::Map1Dto3D>(new AMP::Operator::Map1Dto3D( mapflowcladParams ));
        map1DFlowTo3DFlow2 = boost::shared_ptr<AMP::Operator::Map1Dto3D>(new AMP::Operator::Map1Dto3D( mapflowcladParams ));

        mapCladTo1DFlow1->setZLocations( map1DFlowTo3DFlow1->getZLocations());
        mapCladTo1DFlow2->setZLocations( map1DFlowTo3DFlow2->getZLocations());
        size_t flowVecSize = map1DFlowTo3DFlow1->getNumZlocations();

        boost::shared_ptr<AMP::InputDatabase> mapFlowToClad_db = boost::dynamic_pointer_cast<AMP::InputDatabase>(input_db->getDatabase("MapFlowtoClad"));
        map3DFlowToClad = boost::dynamic_pointer_cast<AMP::Operator::MapSurface>( 
            AMP::Operator::OperatorBuilder::createOperator(meshAdapter2, meshAdapter2, meshAdapter2->getComm(), mapFlowToClad_db ) );
        map3DFlowToClad->setVector(thermalMapToCladVec);
        //------------------------------------------

        cladVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , mapCladTo1DFlow1->getOutputVariable() );
        flowSol1DVec = AMP::LinearAlgebra::SimpleVector::create( flowVecSize , mapCladTo1DFlow1->getOutputVariable() );

        mapCladTo1DFlow1->setVector(cladVec); 
        mapCladTo1DFlow2->setVector(cladVec); 

        AMP::LinearAlgebra::Variable::shared_ptr flowVariable = map1DFlowTo3DFlow1->getOutputVariable() ;

        flowSolVec = AMP::LinearAlgebra::createVector( flowNodalScalarDOF, flowVariable, true );
        flowRhsVec = AMP::LinearAlgebra::createVector( flowNodalScalarDOF, flowVariable, true );
        flowResVec = AMP::LinearAlgebra::createVector( flowNodalScalarDOF, flowVariable, true );

        flowOperator->setZLocations( map1DFlowTo3DFlow1->getZLocations());
        flowJacobian->setZLocations( map1DFlowTo3DFlow1->getZLocations());

        flowOperator->setVector(cladVec);
        flowJacobian->setVector(cladVec);

        flowSolVec->setToScalar(300.0);
        //flowRhsVec->setToScalar(300.0);
        //flowResVec->setToScalar(300.0);
    }


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

    //  AMP::LinearAlgebra::Vector::shared_ptr globalResMultiVectorView = AMP::LinearAlgebra::MultiVector::view( globalResMultiVector, globalComm );
    //---------------------------------------------------------------------------------------------------------------------//

    //------------------------------------------


    //------------------------------------------

    boost::shared_ptr<AMP::InputDatabase> tmp1_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::OperatorParameters> columnMapsParams (new AMP::Operator::OperatorParameters(tmp1_db));
    boost::shared_ptr<AMP::Operator::ColumnOperator>     columnMapstoCladOperator(new AMP::Operator::ColumnOperator(columnMapsParams));
    columnMapstoCladOperator->append(mapPelletAndClad);
    if ( map3DFlowToClad.get() != NULL )
        columnMapstoCladOperator->append(map3DFlowToClad);

    //------------------------------------------
    boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp1;
    boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp2;
    boost::shared_ptr<AMP::Operator::RobinVectorCorrection> robinBoundaryOp3;
    boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters1;
    boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters2;
    boost::shared_ptr<AMP::Operator::NeumannVectorCorrectionParameters> correctionParameters3;
    if ( thermalNonlinearOperator1.get() != NULL ) {
        robinBoundaryOp1 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( thermalNonlinearOperator1->getBoundaryOperator()  );
        correctionParameters1 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp1->getParameters());
        robinBoundaryOp1->setVariableFlux( thermalMapCladToPelletVec );
        robinBoundaryOp1->reset(correctionParameters1);
    }
    if ( thermalNonlinearOperator2.get() != NULL ) {
        robinBoundaryOp2 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( thermalNonlinearOperator2->getBoundaryOperator() ) )->getBoundaryOperator(0) );
        robinBoundaryOp3 = boost::dynamic_pointer_cast<AMP::Operator::RobinVectorCorrection>( (boost::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator> ( thermalNonlinearOperator2->getBoundaryOperator() ) )->getBoundaryOperator(1) );
        correctionParameters2 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp2->getParameters());
        correctionParameters3 = boost::dynamic_pointer_cast<AMP::Operator::NeumannVectorCorrectionParameters>(robinBoundaryOp3->getParameters());
        robinBoundaryOp2->setVariableFlux( thermalMapToCladVec );
        robinBoundaryOp3->setVariableFlux( thermalMapToCladVec );
        robinBoundaryOp2->reset(correctionParameters2);
        robinBoundaryOp3->reset(correctionParameters3);
    }

    //--------------------------------------

    AMP_INSIST(input_db->keyExists("NonlinearSolver"),   "Key ''NonlinearSolver'' is missing!");

    //-------------------------------------
    //Coupling Map to the Nonlinear Operators
    boost::shared_ptr<AMP::InputDatabase> tmp_db (new AMP::InputDatabase("Dummy"));
    boost::shared_ptr<AMP::Operator::OperatorParameters> columnNonlinearParams (new AMP::Operator::OperatorParameters(tmp_db));
    boost::shared_ptr<AMP::Operator::ColumnOperator>     columnNonlinearOperator(new AMP::Operator::ColumnOperator(columnNonlinearParams));
    if ( thermalNonlinearOperator1.get() != NULL ) {
        boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams1(new AMP::Operator::CoupledOperatorParameters(tmp_db));
        coupledNonlinearParams1->d_MapOperator = mapPelletAndClad;
        coupledNonlinearParams1->d_BVPOperator = thermalNonlinearOperator1;
        boost::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearOperator1(new AMP::Operator::CoupledOperator(coupledNonlinearParams1));
        columnNonlinearOperator->append(coupledNonlinearOperator1);
    }
    if ( thermalNonlinearOperator2.get() != NULL ) {
        boost::shared_ptr<AMP::Operator::CoupledOperatorParameters> coupledNonlinearParams2(new AMP::Operator::CoupledOperatorParameters(tmp_db));
        coupledNonlinearParams2->d_MapOperator = columnMapstoCladOperator;
        coupledNonlinearParams2->d_BVPOperator = thermalNonlinearOperator2;
        boost::shared_ptr<AMP::Operator::CoupledOperator> coupledNonlinearOperator2(new AMP::Operator::CoupledOperator(coupledNonlinearParams2));
        columnNonlinearOperator->append(coupledNonlinearOperator2);
    }
    if ( flowOperator.get() != NULL ) {
        boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperatorParameters> coupledNonlinearParams3(new AMP::Operator::CoupledFlowFrapconOperatorParameters(tmp_db));
        coupledNonlinearParams3->d_Map3to1 = mapCladTo1DFlow1 ;
        coupledNonlinearParams3->d_FlowOperator = flowOperator;
        coupledNonlinearParams3->d_Map1to3 = map1DFlowTo3DFlow1;
        coupledNonlinearParams3->d_Mesh = meshAdapter2;
        boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperator>           coupledNonlinearOperator3(new AMP::Operator::CoupledFlowFrapconOperator(coupledNonlinearParams3));
        columnNonlinearOperator->append(coupledNonlinearOperator3);
    }

    //---------------------------------------
    //Column of Coupled Operators
    boost::shared_ptr<AMP::Operator::OperatorParameters> columnLinearParams(new AMP::Operator::OperatorParameters(tmp_db));
    boost::shared_ptr<AMP::Operator::ColumnOperator>     coupledLinearOperator(new AMP::Operator::ColumnOperator(columnLinearParams));
    if ( thermalLinearOperator1.get() != NULL )
        coupledLinearOperator->append(thermalLinearOperator1);
    if ( thermalLinearOperator2.get() != NULL )
        coupledLinearOperator->append(thermalLinearOperator2);
    boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperator> coupledlinearOperator3;
    if ( flowJacobian.get() != NULL ) {
        boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperatorParameters> coupledlinearParams3(new AMP::Operator::CoupledFlowFrapconOperatorParameters(tmp_db));
        coupledlinearParams3->d_Map3to1 = mapCladTo1DFlow2 ;
        coupledlinearParams3->d_FlowOperator = flowJacobian;
        coupledlinearParams3->d_Map1to3 = map1DFlowTo3DFlow2;
        coupledlinearParams3->d_Mesh = meshAdapter2;
        coupledlinearOperator3 = boost::shared_ptr<AMP::Operator::CoupledFlowFrapconOperator>(new AMP::Operator::CoupledFlowFrapconOperator(coupledlinearParams3));
        coupledLinearOperator->append(coupledlinearOperator3);
    }

    //---------------------------------------------------------------------------------------------------------------------//
    if ( neutronicsOperator.get() != NULL ) {
        neutronicsOperator->setTimeStep(0);
        neutronicsOperator->apply(nullVec, nullVec, specificPowerGpVec, 1., 0.);
        thermalRhsVec1->zero();
        specificPowerGpVecToPowerDensityNodalVecOperatator->apply(nullVec, specificPowerGpVec, thermalRhsVec1, 1., 0.);
    }

    //We need to reset the linear operator before the solve since TrilinosML does
    //the factorization of the matrix during construction and so the matrix must
    //be correct before constructing the TrilinosML object.
    //The thermal operator does not expect an apply to be called before calling
    //getJacobianParams and so it need not be called. So, any of the following
    //apply calls will work:
    coupledLinearOperator->reset(columnNonlinearOperator->getJacobianParameters(globalSolMultiVector));
    columnNonlinearOperator->apply(nullVec, globalSolMultiVector, globalResMultiVector, 1.0, 0.0);
    AMP::pout<<"Initial Global Residual Norm: "<<std::setprecision(12)<<globalResMultiVector->L2Norm()<<std::endl;
    AMP::pout<<"Initial Temperature Residual Norm: "<<std::setprecision(12)<<globalResVec->L2Norm()<<std::endl;
    if ( flowResVec.get() != NULL )
        AMP::pout<<"Initial Flow Residual Norm: "<<std::setprecision(12)<<flowResVec->L2Norm()<<std::endl;

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
    if ( thermalLinearOperator1.get() != NULL ) {
        boost::shared_ptr<AMP::Database> thermalPreconditioner_db1 = columnPreconditioner_db->getDatabase("pelletThermalPreconditioner");
        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams1(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db1));
        thermalPreconditionerParams1->d_pOperator = thermalLinearOperator1;
        boost::shared_ptr<AMP::Solver::TrilinosMLSolver> thermalPreconditioner1(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams1));
        columnPreconditioner->append(thermalPreconditioner1);
    }
    if ( thermalLinearOperator2.get() != NULL ) {
        boost::shared_ptr<AMP::Database> thermalPreconditioner_db2 = columnPreconditioner_db->getDatabase("cladThermalPreconditioner");
        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> thermalPreconditionerParams2(new AMP::Solver::SolverStrategyParameters(thermalPreconditioner_db2));
        thermalPreconditionerParams2->d_pOperator = thermalLinearOperator2;
        boost::shared_ptr<AMP::Solver::TrilinosMLSolver> thermalPreconditioner2(new AMP::Solver::TrilinosMLSolver(thermalPreconditionerParams2));
        columnPreconditioner->append(thermalPreconditioner2);
    }
    if ( flowJacobian.get() != NULL ) {
        boost::shared_ptr<AMP::Database> JacobianSolver_db = input_db->getDatabase("Flow1DSolver"); 
        boost::shared_ptr<AMP::Solver::SolverStrategyParameters> flowSolverParams (new AMP::Solver::SolverStrategyParameters(JacobianSolver_db));
        flowSolverParams->d_pOperator = flowJacobian;
        boost::shared_ptr<AMP::Solver::Flow1DSolver> flowJacobianSolver(new AMP::Solver::Flow1DSolver(flowSolverParams));

        boost::shared_ptr<AMP::InputDatabase> CoupledJacobianSolver_db (new AMP::InputDatabase("Dummy"));
        CoupledJacobianSolver_db->putInteger("max_iterations",1); 
        CoupledJacobianSolver_db->putDouble("max_error",1.0e-6); 
        boost::shared_ptr<AMP::Solver::CoupledFlow1DSolverParameters> coupledFlowSolverParams (new AMP::Solver::CoupledFlow1DSolverParameters(CoupledJacobianSolver_db ));
        coupledFlowSolverParams->d_flow1DSolver = flowJacobianSolver;
        coupledFlowSolverParams->d_pOperator = boost::dynamic_pointer_cast< AMP::Operator::Operator>(coupledlinearOperator3);
        boost::shared_ptr<AMP::Solver::CoupledFlow1DSolver> CoupledFlowJacobianSolver(new AMP::Solver::CoupledFlow1DSolver(coupledFlowSolverParams));
        columnPreconditioner->append(CoupledFlowJacobianSolver);
    }

    //--------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    boost::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver = nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner(columnPreconditioner);

    //-------------------------------------
    nonlinearSolver->setZeroInitialGuess(false);

#ifdef USE_EXT_SILO
    // Register the quantities to plot
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter("Silo");
    if ( meshAdapter1.get()!=NULL ) {
        siloWriter->registerVector( globalSolVec, meshAdapter1, AMP::Mesh::Vertex, "PelletTemperature" );
    }
    if ( meshAdapter2.get()!=NULL ) {
        siloWriter->registerVector( globalSolVec, meshAdapter2, AMP::Mesh::Vertex, "CladTemperature" );
        siloWriter->registerVector( flowSolVec, surfaceMesh, AMP::Mesh::Vertex, "FlowTemperature" );
    }
#endif

    for ( int tstep = 0; tstep < 1; tstep++ ) {
        if ( neutronicsOperator.get()!=NULL ) {
            neutronicsOperator->setTimeStep(tstep);
            neutronicsOperator->apply(nullVec, nullVec, specificPowerGpVec, 1., 0.);
        }
        
        if ( robinBoundaryOp1.get() != NULL )
            robinBoundaryOp1->reset(correctionParameters1);
        if ( robinBoundaryOp2.get() != NULL )
            robinBoundaryOp2->reset(correctionParameters2);
        if ( robinBoundaryOp3.get() != NULL )
            robinBoundaryOp3->reset(correctionParameters3);

        if ( meshAdapter1.get()!=NULL ) {
            thermalRhsVec1->zero();
            // specificPowerGpVec is in Watts/kilogram
            specificPowerGpVecToPowerDensityNodalVecOperatator->apply(nullVec, specificPowerGpVec, thermalRhsVec1, 1., 0.);
        }
        if ( thermalNonlinearOperator1.get()!=NULL ) {
            thermalNonlinearOperator1->modifyRHSvector(thermalRhsVec1);
            thermalNonlinearOperator1->modifyInitialSolutionVector(thermalSolVec1);
        }
        if ( thermalNonlinearOperator2.get()!=NULL ) {
            thermalNonlinearOperator2->modifyRHSvector(thermalRhsVec2);
            thermalNonlinearOperator2->modifyInitialSolutionVector(thermalSolVec2);
        }

        AMP::pout<<"Initial Guess  Norm for Step " << tstep << " is: "<<globalSolMultiVector->L2Norm()<<std::endl;
        AMP::pout<<"Initial Guess  Norm12 for Step " << tstep << " is: "<<globalSolVec->L2Norm()<<std::endl;
        if ( flowSolVec.get() != NULL )
            AMP::pout<<"Initial Guess  Flow   for Step " << tstep << " is: "<<flowSolVec->L2Norm()<<std::endl;
        AMP::pout<<"Initial Source Norm for Step " << tstep << " is: "<<globalRhsMultiVector->L2Norm()<<std::endl;
        if ( thermalRhsVec1.get() != NULL ) {
            AMP::pout<<"Initial Source Norm1 for Step " << tstep << " is: "<<thermalRhsVec1->L2Norm()<<std::endl;
            AMP::pout<<"Initial Guess  Norm1 for Step " << tstep << " is: "<<thermalSolVec1->L2Norm()<<std::endl;
            AMP::pout<<"Initial Power  Norm1 for Step " << tstep << " is: "<<specificPowerGpVec->L2Norm()<<std::endl;
        }
        if ( thermalRhsVec2.get() != NULL ) {
            AMP::pout<<"Initial Source Norm2 for Step " << tstep << " is: "<<thermalRhsVec2->L2Norm()<<std::endl;
            AMP::pout<<"Initial Guess  Norm2 for Step " << tstep << " is: "<<thermalSolVec2->L2Norm()<<std::endl;
        }
        globalResMultiVector->zero();
        columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
        AMP::pout<<"Initial Global Residual Norm for Step " << tstep << " is: "<<globalResMultiVector->L2Norm()<<std::endl;
        AMP::pout<<"Initial Temperature Residual Norm for Step " << tstep << " is: "<<globalResVec->L2Norm()<<std::endl;
        if ( flowResVec.get() != NULL )
            AMP::pout<<"Initial Flow Residual for Step " << tstep << " is: "<<flowResVec->L2Norm()<<std::endl;

        nonlinearSolver->solve(globalRhsMultiVectorView, globalSolMultiVectorView);

        columnNonlinearOperator->apply(globalRhsMultiVector, globalSolMultiVector, globalResMultiVector, 1.0, -1.0);
        AMP::pout<<"Final   Residual Norm for Step " << tstep << " is: "<<globalResMultiVector->L2Norm()<<std::endl;
        #ifdef USE_EXT_SILO
            siloWriter->writeFile( silo_name , tstep );
        #endif

        if ( thermalSolVec1.get() != NULL )
            std::cout<< "The Fuel Max value - " << thermalSolVec1->max() << std::endl;
        if ( thermalSolVec2.get() != NULL )
            std::cout<< "The Clad Min value - " << thermalSolVec2->min() << std::endl;
        if ( flowSolVec.get() != NULL )
            std::cout << "Flow Max:Min values -  " << flowSolVec->max() << " " << flowSolVec->min() << std::endl;

        if ( meshAdapter2.get()!=NULL ) {
            std::cout<<"Intermediate Flow Solution " <<std::endl;
            mapCladTo1DFlow1->setVector(flowSol1DVec); 
	        mapCladTo1DFlow1->apply(nullVec, thermalMapToCladVec , nullVec, 1, 0);
            size_t flowVecSize = map1DFlowTo3DFlow1->getNumZlocations();
            std::vector<double> expectedSolution(flowVecSize, 0);
            expectedSolution = (input_db->getDatabase("regression"))->getDoubleArray("expectedSolution");
            for(unsigned int i=0; i<flowVecSize ; i++) {
                if( !AMP::Utilities::approx_equal(expectedSolution[i], flowSol1DVec->getValueByLocalID(i), 1e-6)) {
                    if ( meshAdapter2->getComm().getRank() == 0 ){
                        printf("solution: %.7e expected: %.7e \n", 
                        flowSol1DVec->getValueByLocalID(i),expectedSolution[i] );
                    
                    }
                    ut->failure("solution is different for "+silo_name);
                }
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


