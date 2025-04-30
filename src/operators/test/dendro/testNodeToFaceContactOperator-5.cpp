#include "AMP/IO/PIO.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/Mesh.h"
#include "AMP/mesh/MeshFactory.h"
#include "AMP/mesh/euclidean_geometry_tools.h"
#include "AMP/mesh/latex_visualization_tools.h"
#include "AMP/mesh/libmesh/ReadTestMesh.h"
#include "AMP/mesh/libmesh/libmeshMesh.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/TrilinosMatrixShellOperator.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/contact/NodeToGeomType::FaceContactOperator.h"
#include "AMP/operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "AMP/operators/mechanics/IsotropicElasticModel.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsMaterialModel.h"
#include "AMP/operators/mechanics/MechanicsModelParameters.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "externVars.h"
#include "testNodeToGeomType::FaceContactOperator.h"

#include <fstream>
#include <set>


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    int rank = globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_driver_" + std::to_string( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    // Load the input file
    globalComm.barrier();
    double inpReadBeginTime = MPI_Wtime();


    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    globalComm.barrier();
    double inpReadEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished parsing the input file in " << ( inpReadEndTime - inpReadBeginTime )
                  << " seconds." << std::endl;
    }

    bool useML = input_db->getWithDefault<bool>( "useML", false );
    bool cladExpansionConstrained =
        input_db->getWithDefault<bool>( "cladExpansionConstrained", true );
    bool useLevitatingFuel = input_db->getWithDefault<bool>( "useLevitatingFuel", true );
    std::string prefixFileName =
        input_db->getWithDefault<std::string>( "prefixFileName", "TATA_0" );
    double scaleSolution       = input_db->getWithDefault<double>( "scaleSolution", 1.0 );
    double cladNeedALittleHelp = input_db->getWithDefault<double>( "cladNeedALittleHelp", 0.0 );
    double fuelNeedALittleHelp = input_db->getWithDefault<double>( "fuelNeedALittleHelp", -1.0 );
    bool contactIsFrictionless = input_db->getWithDefault<bool>( "contactIsFrictionless", false );
    double shrinkFactor        = input_db->getWithDefault<double>( "shrinkFactor", 0.0 );

    // Load the meshes
    globalComm.barrier();
    double meshBeginTime = MPI_Wtime();

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "Mesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( globalComm );
    auto meshAdapter = AMP::Mesh::MeshFactory::create( meshParams );

    globalComm.barrier();
    double meshEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished reading the mesh in " << ( meshEndTime - meshBeginTime )
                  << " seconds." << std::endl;
    }

    // Create a DOF manager
    int dofsPerNode     = 3;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto dispDofManager = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, dofsPerNode, split );

    // Build a column operator and a column preconditioner
    auto columnOperator          = std::make_shared<AMP::Operator::ColumnOperator>();
    auto linearSolver_db         = input_db->getDatabase( "LinearSolver" );
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        std::make_shared<AMP::Solver::ColumnSolverParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = columnOperator;
    auto columnPreconditioner =
        std::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    // Get the mechanics material models for the contact operator and for computing stresses
    auto fuelModel_db = input_db->getDatabase( "FuelMechanicsMaterialModel" );
    auto fuelMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( fuelModel_db );
    auto fuelMechanicsMaterialModel =
        std::make_shared<AMP::Operator::IsotropicElasticModel>( fuelMechanicsMaterialModelParams );
    auto cladModel_db = input_db->getDatabase( "CladMechanicsMaterialModel" );
    auto cladMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( cladModel_db );
    auto cladMechanicsMaterialModel =
        std::make_shared<AMP::Operator::IsotropicElasticModel>( cladMechanicsMaterialModelParams );

    // Build the contact operators
    auto bottomPelletTopPelletContact_db =
        input_db->getDatabase( "BottomPelletTopPelletContactOperator" );
    auto bottomPelletTopPelletContactOperatorParams =
        std::make_shared<AMP::Operator::ContactOperatorParameters>(
            bottomPelletTopPelletContact_db );
    bottomPelletTopPelletContactOperatorParams->d_DOFsPerNode = dofsPerNode;
    bottomPelletTopPelletContactOperatorParams->d_DOFManager  = dispDofManager;
    bottomPelletTopPelletContactOperatorParams->d_GlobalComm  = globalComm;
    bottomPelletTopPelletContactOperatorParams->d_Mesh        = meshAdapter;
    bottomPelletTopPelletContactOperatorParams->d_MasterMechanicsMaterialModel =
        fuelMechanicsMaterialModel;
    bottomPelletTopPelletContactOperatorParams->reset();
    auto bottomPelletTopPelletContactOperator =
        std::make_shared<AMP::Operator::NodeToGeomType::FaceContactOperator>(
            bottomPelletTopPelletContactOperatorParams );
    bottomPelletTopPelletContactOperator->initialize();
    bottomPelletTopPelletContactOperator->setContactIsFrictionless( contactIsFrictionless );
    // bottomPelletTopPelletContactOperator->setContactIsFrictionless(bottomPelletTopPelletContact_db->getWithDefault<bool>("ContactIsFrictionless",false));

    auto bottomPelletCladContact_db = input_db->getDatabase( "BottomPelletCladContactOperator" );
    auto bottomPelletCladContactOperatorParams =
        std::make_shared<AMP::Operator::ContactOperatorParameters>( bottomPelletCladContact_db );
    bottomPelletCladContactOperatorParams->d_DOFsPerNode = dofsPerNode;
    bottomPelletCladContactOperatorParams->d_DOFManager  = dispDofManager;
    bottomPelletCladContactOperatorParams->d_GlobalComm  = globalComm;
    bottomPelletCladContactOperatorParams->d_Mesh        = meshAdapter;
    bottomPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel =
        fuelMechanicsMaterialModel;
    bottomPelletCladContactOperatorParams->reset();
    auto bottomPelletCladContactOperator =
        std::make_shared<AMP::Operator::NodeToGeomType::FaceContactOperator>(
            bottomPelletCladContactOperatorParams );
    bottomPelletCladContactOperator->initialize();
    bottomPelletCladContactOperator->setContactIsFrictionless( contactIsFrictionless );
    // bottomPelletCladContactOperator->setContactIsFrictionless(bottomPelletCladContact_db->getWithDefault<bool>("ContactIsFrictionless",false));

    auto topPelletCladContact_db = input_db->getDatabase( "TopPelletCladContactOperator" );
    auto topPelletCladContactOperatorParams =
        std::make_shared<AMP::Operator::ContactOperatorParameters>( topPelletCladContact_db );
    topPelletCladContactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    topPelletCladContactOperatorParams->d_DOFManager                   = dispDofManager;
    topPelletCladContactOperatorParams->d_GlobalComm                   = globalComm;
    topPelletCladContactOperatorParams->d_Mesh                         = meshAdapter;
    topPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel = fuelMechanicsMaterialModel;
    topPelletCladContactOperatorParams->reset();
    auto topPelletCladContactOperator =
        std::make_shared<AMP::Operator::NodeToGeomType::FaceContactOperator>(
            topPelletCladContactOperatorParams );
    topPelletCladContactOperator->initialize();
    topPelletCladContactOperator->setContactIsFrictionless( contactIsFrictionless );
    // topPelletCladContactOperator->setContactIsFrictionless(topPelletCladContact_db->getWithDefault<bool>("ContactIsFrictionless",false));

    // Build the BVP operators
    std::shared_ptr<AMP::Operator::LinearBVPOperator> bottomPelletBVPOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> topPelletBVPOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> cladBVPOperator;


    AMP::Mesh::MeshID bottomPelletMeshID = bottomPelletTopPelletContactOperator->getMasterMeshID();
    AMP_ASSERT( bottomPelletMeshID == bottomPelletCladContactOperator->getMasterMeshID() );
    auto bottomPelletMeshAdapter = meshAdapter->Subset( bottomPelletMeshID );
    if ( bottomPelletMeshAdapter ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> bottomPelletElementPhysicsModel;
        bottomPelletBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( bottomPelletMeshAdapter,
                                                            "BottomPelletBVPOperator",
                                                            input_db,
                                                            bottomPelletElementPhysicsModel ) );
        columnOperator->append( bottomPelletBVPOperator );

        if ( !useML ) {
            auto bottomPelletSolver_db = columnPreconditioner_db->getDatabase( "DummySolver" );
            auto bottomPelletSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( bottomPelletSolver_db );
            bottomPelletSolverParams->d_pOperator = bottomPelletBVPOperator;
            bottomPelletSolverParams->d_comm      = bottomPelletMeshAdapter->getComm();
            auto bottomPelletSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( bottomPelletSolverParams );
            columnPreconditioner->append( bottomPelletSolver );
        } else {
            auto bottomPelletSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto bottomPelletSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( bottomPelletSolver_db );
            bottomPelletSolverParams->d_pOperator = bottomPelletBVPOperator;
            auto bottomPelletSolver =
                std::make_shared<AMP::Solver::TrilinosMLSolver>( bottomPelletSolverParams );
            columnPreconditioner->append( bottomPelletSolver );
        } // end if
    }     // end if

    auto topPelletMeshID = bottomPelletTopPelletContactOperator->getSlaveMeshID();
    AMP_ASSERT( topPelletMeshID == topPelletCladContactOperator->getMasterMeshID() );
    auto topPelletMeshAdapter = meshAdapter->Subset( topPelletMeshID );
    if ( topPelletMeshAdapter ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> topPelletElementPhysicsModel;
        topPelletBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( topPelletMeshAdapter,
                                                            "TopPelletBVPOperator",
                                                            input_db,
                                                            topPelletElementPhysicsModel ) );
        columnOperator->append( topPelletBVPOperator );

        if ( !useML ) {
            auto topPelletSolver_db = columnPreconditioner_db->getDatabase( "DummySolver" );
            auto topPelletSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( topPelletSolver_db );
            topPelletSolverParams->d_pOperator = topPelletBVPOperator;
            topPelletSolverParams->d_comm      = topPelletMeshAdapter->getComm();
            auto topPelletSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( topPelletSolverParams );
            columnPreconditioner->append( topPelletSolver );
        } else {
            auto topPelletSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto topPelletSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( topPelletSolver_db );
            topPelletSolverParams->d_pOperator = topPelletBVPOperator;
            auto topPelletSolver =
                std::make_shared<AMP::Solver::TrilinosMLSolver>( topPelletSolverParams );
            columnPreconditioner->append( topPelletSolver );
        } // end if
    }     // end if

    auto cladMeshID = bottomPelletCladContactOperator->getSlaveMeshID();
    AMP_ASSERT( cladMeshID == topPelletCladContactOperator->getSlaveMeshID() );
    auto cladMeshAdapter = meshAdapter->Subset( cladMeshID );
    if ( cladMeshAdapter ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> cladElementPhysicsModel;
        cladBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                cladMeshAdapter, "CladBVPOperator", input_db, cladElementPhysicsModel ) );
        columnOperator->append( cladBVPOperator );

        if ( useML ) {
            auto cladSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto cladSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( cladSolver_db );
            cladSolverParams->d_pOperator = cladBVPOperator;
            auto cladSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( cladSolverParams );
            columnPreconditioner->append( cladSolver );
        } else {
            auto cladSolver_db = columnPreconditioner_db->getDatabase( "DummySolver" );
            auto cladSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( cladSolver_db );
            cladSolverParams->d_pOperator = cladBVPOperator;
            cladSolverParams->d_comm      = cladMeshAdapter->getComm();
            auto cladSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( cladSolverParams );
            columnPreconditioner->append( cladSolver );
        } // end if

    } // end if


    {
        auto bottomPelletTopPelletContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        auto bottomPelletTopPelletContactPreconditionerParams =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
                bottomPelletTopPelletContactPreconditioner_db );
        bottomPelletTopPelletContactPreconditionerParams->d_pOperator =
            bottomPelletTopPelletContactOperator;
        auto bottomPelletTopPelletContactPreconditioner =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolver>(
                bottomPelletTopPelletContactPreconditionerParams );
        columnPreconditioner->append( bottomPelletTopPelletContactPreconditioner );
    }
    {
        auto bottomPelletCladContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        auto bottomPelletCladContactPreconditionerParams =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
                bottomPelletCladContactPreconditioner_db );
        bottomPelletCladContactPreconditionerParams->d_pOperator = bottomPelletCladContactOperator;
        auto bottomPelletCladContactPreconditioner =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolver>(
                bottomPelletCladContactPreconditionerParams );
        columnPreconditioner->append( bottomPelletCladContactPreconditioner );
    }
    {
        auto topPelletCladContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        auto topPelletCladContactPreconditionerParams =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
                topPelletCladContactPreconditioner_db );
        topPelletCladContactPreconditionerParams->d_pOperator = topPelletCladContactOperator;
        auto topPelletCladContactPreconditioner =
            std::make_shared<AMP::Solver::ConstraintsEliminationSolver>(
                topPelletCladContactPreconditionerParams );
        columnPreconditioner->append( topPelletCladContactPreconditioner );
    }

    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> bottomPelletConstraints;
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> topPelletConstraints;
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> cladConstraints;
    if ( !useLevitatingFuel ) {
        double fuelOuterRadius = input_db->getScalar<double>( "FuelOuterRadius" );
        double dishRadius      = input_db->getWithDefault<double>( "DishRadius", 0.0 );
        makeConstraintsOnFuel(
            bottomPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
            fuelOuterRadius,
            bottomPelletConstraints,
            true,
            dishRadius );
        makeConstraintsOnFuel(
            bottomPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
            fuelOuterRadius,
            bottomPelletConstraints,
            false );
        makeConstraintsOnFuel(
            topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
            fuelOuterRadius,
            topPelletConstraints,
            true,
            dishRadius );
        makeConstraintsOnFuel(
            topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
            fuelOuterRadius,
            topPelletConstraints,
            false );
    } // end if
    if ( !cladExpansionConstrained ) {
        double cladInnerRadius = input_db->getScalar<double>( "CladInnerRadius" );
        double cladOuterRadius = input_db->getScalar<double>( "CladOuterRadius" );
        double cladHeight      = input_db->getScalar<double>( "CladHeight" );
        makeConstraintsOnClad( cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex ),
                               cladInnerRadius,
                               cladOuterRadius,
                               cladHeight,
                               cladConstraints );
    } // end if

    // Items for computing the RHS correction due to thermal expansion
    auto fuelTemperatureRhs_db = input_db->getDatabase( "FuelTemperatureRHSVectorCorrection" );
    auto cladTemperatureRhs_db = input_db->getDatabase( "CladTemperatureRHSVectorCorrection" );
    auto tempVar               = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );
    auto dispVar               = columnOperator->getOutputVariable();
    auto tempDofManager        = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, 1, split );
    auto tempVec    = AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    auto refTempVec = tempVec->clone();
    auto sigma_xx   = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xx" ), split );
    auto sigma_yy = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_yy" ), split );
    auto sigma_zz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_zz" ), split );
    auto sigma_yz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_yz" ), split );
    auto sigma_xz = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xz" ), split );
    auto sigma_xy = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_xy" ), split );
    auto sigma_eff = AMP::LinearAlgebra::createVector(
        tempDofManager, std::make_shared<AMP::LinearAlgebra::Variable>( "sigma_eff" ), split );

    globalComm.barrier();
    double tempCompBeginTime = MPI_Wtime();

    double dummyFuelThermalConductivity = 0.0; // not used k_f = a+b/(c+Td)
    double linearHeatGenerationRate     = input_db->getScalar<double>( "LinearHeatGenerationRate" );
    double fuelOuterRadius              = input_db->getScalar<double>( "FuelOuterRadius" );
    double cladInnerRadius              = input_db->getScalar<double>( "CladInnerRadius" );
    double cladOuterRadius              = input_db->getScalar<double>( "CladOuterRadius" );
    double cladThermalConductivity      = input_db->getScalar<double>( "CladThermalConductivity" );

    double gapThermalConductivity = input_db->getScalar<double>( "GapThermalConductivity" );
    double moderatorTemperature   = input_db->getScalar<double>( "ModeratorTemperature" );
    double moderatorHeatTransferCoefficient =
        input_db->getScalar<double>( "ModeratorHeatTransferCoefficient" );

    double referenceTemperature = input_db->getScalar<double>( "ReferenceTemperature" );

    double cladOuterRadiusTemperature =
        moderatorTemperature + linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                   ( cladOuterRadius * moderatorHeatTransferCoefficient );
    double cladInnerRadiusTemperature =
        cladOuterRadiusTemperature + linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                         cladThermalConductivity *
                                         std::log( cladOuterRadius / cladInnerRadius );
    double gapAverageRadius = 0.5 * ( cladInnerRadius + fuelOuterRadius );
    double gapHeatTranferCoefficient =
        gapThermalConductivity / ( cladInnerRadius - fuelOuterRadius );
    double fuelOuterRadiusTemperature =
        cladInnerRadiusTemperature + linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                         ( gapAverageRadius * gapHeatTranferCoefficient );
    double fuelCenterLineTemperature = fuelOuterRadiusTemperature;
    newton_solver_t<double> solver;
    solver.set( &my_f, &my_ijmf );
    my_p.r       = 0.0;
    my_p.q_prime = linearHeatGenerationRate;
    my_p.T_fo    = fuelOuterRadiusTemperature;
    my_p.R_fo    = fuelOuterRadius;
    solver.apply( fuelCenterLineTemperature, &my_p );
    if ( !rank ) {
        std::cout << "cladOuterRadiusTemperature=" << cladOuterRadiusTemperature << "\n";
        std::cout << "cladInnerRadiusTemperature=" << cladInnerRadiusTemperature << "\n";
        std::cout << "fuelOuterRadiusTemperature=" << fuelOuterRadiusTemperature << "\n";
        std::cout << "fuelCenterLineTemperature=" << fuelCenterLineTemperature << "\n";
    }

    refTempVec->setToScalar( referenceTemperature );
    tempVec->setToScalar( referenceTemperature );

    computeFuelTemperature( bottomPelletMeshAdapter,
                            tempVec,
                            fuelOuterRadius,
                            fuelOuterRadiusTemperature,
                            linearHeatGenerationRate,
                            dummyFuelThermalConductivity );
    computeFuelTemperature( topPelletMeshAdapter,
                            tempVec,
                            fuelOuterRadius,
                            fuelOuterRadiusTemperature,
                            linearHeatGenerationRate,
                            dummyFuelThermalConductivity );
    computeCladTemperature( cladMeshAdapter,
                            tempVec,
                            cladInnerRadius,
                            cladOuterRadius,
                            linearHeatGenerationRate,
                            cladOuterRadiusTemperature,
                            cladThermalConductivity );

    globalComm.barrier();
    double tempCompEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished computing the temperature profile in "
                  << ( tempCompEndTime - tempCompBeginTime ) << " seconds." << std::endl;
    }

    double fuelThermalExpansionCoefficient =
        ( fuelTemperatureRhs_db->getDatabase( "RhsMaterialModel" ) )
            ->getScalar<double>( "THERMAL_EXPANSION_COEFFICIENT" );
    double cladThermalExpansionCoefficient =
        ( cladTemperatureRhs_db->getDatabase( "RhsMaterialModel" ) )
            ->getScalar<double>( "THERMAL_EXPANSION_COEFFICIENT" );
    bottomPelletTopPelletContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );
    bottomPelletCladContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );
    topPelletCladContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnVar    = columnOperator->getOutputVariable();
    auto columnSolVec = AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    auto columnRhsVec = AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();
    AMP::LinearAlgebra::Vector::shared_ptr bottomPelletCor;
    AMP::LinearAlgebra::Vector::shared_ptr topPelletCor;
    AMP::LinearAlgebra::Vector::shared_ptr cladCor;

    auto activeSetBeforeUpdateVec = sigma_eff->clone();
    auto activeSetAfterUpdateVec  = sigma_eff->clone();
    auto contactPressureVec       = sigma_eff->clone();
    auto surfaceTractionVec       = columnSolVec->clone();
    auto normalVectorVec          = columnSolVec->clone();

    if ( shrinkFactor != 0.0 ) {
        AMP_ASSERT( ( shrinkFactor > 0.0 ) && ( shrinkFactor < 1.0 ) );
        shrinkMesh( cladMeshAdapter, shrinkFactor );
    }

    if ( fuelNeedALittleHelp != -1.0 ) {
        {
            auto it = topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 );
            auto it_begin = it.begin();
            auto it_end   = it.end();
            std::vector<size_t> dofs;
            double epsilon = 1.0e-12;
            double radius;
            for ( it = it_begin; it != it_end; ++it ) {
                auto coord = it->coord();
                radius     = std::sqrt( std::pow( coord[0], 2 ) + std::pow( coord[1], 2 ) );
                if ( radius > fuelNeedALittleHelp ) {
                    dispDofManager->getDOFs( it->globalID(), dofs );
                    columnSolVec->setValueByGlobalID( dofs[2], 0.0001 );
                } // end if
            }     // end for
        }
        bottomPelletTopPelletContactOperator->updateActiveSetWithALittleHelp( columnSolVec );
        columnSolVec->zero();
    } else {
        bool skipDisplaceMesh = true;
        bottomPelletTopPelletContactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
    } // end if
    if ( cladNeedALittleHelp != 0.0 ) {
        {
            auto it       = cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
            auto it_begin = it.begin();
            auto it_end   = it.end();
            std::vector<size_t> dofs;
            double epsilon = 1.0e-12;
            for ( it = it_begin; it != it_end; ++it ) {
                dispDofManager->getDOFs( it->globalID(), dofs );
                columnSolVec->setValueByGlobalID( dofs[2], -cladNeedALittleHelp );
            } // end for
        }
        topPelletCladContactOperator->updateActiveSetWithALittleHelp( columnSolVec );
        columnSolVec->zero();
        {
            auto it       = cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
            auto it_begin = it.begin();
            autoit_end    = it.end();
            std::vector<size_t> dofs;
            double epsilon = 1.0e-12;
            for ( it = it_begin; it != it_end; ++it ) {
                dispDofManager->getDOFs( it->globalID(), dofs );
                columnSolVec->setValueByGlobalID( dofs[2], cladNeedALittleHelp );
            } // end for
        }
        bottomPelletCladContactOperator->updateActiveSetWithALittleHelp( columnSolVec );
        columnSolVec->zero();
    } else {
        bool skipDisplaceMesh = true;
        bottomPelletCladContactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
        topPelletCladContactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
    } // end if

    auto contactShiftVec = createVector( dispDofManager, columnVar, split );
    contactShiftVec->zero();
    columnSolVec->zero();
    columnOperator->append( bottomPelletTopPelletContactOperator );
    columnOperator->append( bottomPelletCladContactOperator );
    columnOperator->append( topPelletCladContactOperator );


    // Build a matrix shell operator to use the column operator with the petsc krylov solvers
    auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
    auto matrixShellParams =
        std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase );
    auto matrixShellOperator =
        std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );

    int numBottomPelletLocalNodes = 0;
    int numTopPelletLocalNodes    = 0;
    int numCladLocalNodes         = 0;
    if ( bottomPelletMeshAdapter ) {
        numBottomPelletLocalNodes =
            bottomPelletMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( topPelletMeshAdapter ) {
        numTopPelletLocalNodes =
            topPelletMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( cladMeshAdapter ) {
        numCladLocalNodes = cladMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    int matLocalSize =
        dofsPerNode * ( numBottomPelletLocalNodes + numTopPelletLocalNodes + numCladLocalNodes );
    AMP_ASSERT( matLocalSize == static_cast<int>( dispDofManager->numLocalDOF() ) );
    matrixShellOperator->setComm( globalComm );
    matrixShellOperator->setMatLocalRowSize( matLocalSize );
    matrixShellOperator->setMatLocalColumnSize( matLocalSize );
    matrixShellOperator->setOperator( columnOperator );

    auto linearSolverParams =
        std::make_shared < AMP::Solver::SolverStrategyParameters ? ( linearSolver_db );
    linearSolverParams->d_pOperator     = matrixShellOperator;
    linearSolverParams->d_comm          = globalComm;
    linearSolverParams->d_pNestedSolver = columnPreconditioner;
    auto linearSolver = std::make_shared < AMP::Solver::PetscKrylovSolver ? ( linearSolverParams );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    auto fullThermalLoadingTempMinusRefTempVec = tempVec->clone();
    fullThermalLoadingTempMinusRefTempVec->subtract( tempVec, refTempVec );

    size_t maxActiveSetIterations = input_db->getWithDefault<size_t>( "maxActiveSetIterations", 5 );
    size_t maxThermalLoadingIterations =
        input_db->getWithDefault<size_t>( "maxThermalLoadingIterations", 5 );
    std::vector<double> scalingFactors;
    std::vector<int> maxIterations;
    bool customLoading = input_db->getWithDefault<bool>( "customLoading", false );
    if ( customLoading ) {
        scalingFactors = input_db->getVector<double>( "scalingFactors" );
        maxIterations  = input_db->getVector<int>( "maxIterations" );
        AMP_ASSERT( scalingFactors.size() == maxIterations.size() );
        maxThermalLoadingIterations = scalingFactors.size();
    } // end if

    int TOTO_count = 0;
    for ( size_t thermalLoadingIteration = 0; thermalLoadingIteration < maxThermalLoadingIterations;
          ++thermalLoadingIteration ) {
        double scalingFactor = static_cast<double>( thermalLoadingIteration + 1 ) /
                               static_cast<double>( maxThermalLoadingIterations );
        if ( customLoading ) {
            scalingFactor          = scalingFactors[thermalLoadingIteration];
            maxActiveSetIterations = maxIterations[thermalLoadingIteration];
        } // end if

        if ( !rank ) {
            std::cout << "THERMAL LOADING " << thermalLoadingIteration + 1 << "/"
                      << maxThermalLoadingIterations << "  (" << scalingFactor << ")\n";
        }
        tempVec->axpy( scalingFactor, fullThermalLoadingTempMinusRefTempVec, refTempVec );

        for ( size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations;
              ++activeSetIteration ) {
            if ( !rank ) {
                std::cout << "ACTIVE SET ITERATION #" << activeSetIteration + 1 << std::endl;
            }
            ++TOTO_count;

            columnSolVec->zero();
            columnRhsVec->zero();

            // compute thermal load f
            {
                AMP::LinearAlgebra::VS_Mesh bottomPelletVectorSelector( bottomPelletMeshAdapter );
                auto bottomPelletRhsVec = columnRhsVec->select( bottomPelletVectorSelector );
                computeTemperatureRhsVector( bottomPelletMeshAdapter,
                                             fuelTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             bottomPelletRhsVec );

                AMP::LinearAlgebra::VS_Mesh topPelletVectorSelector( topPelletMeshAdapter );
                auto topPelletRhsVec = columnRhsVec->select( topPelletVectorSelector );
                computeTemperatureRhsVector( topPelletMeshAdapter,
                                             fuelTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             topPelletRhsVec );

                AMP::LinearAlgebra::VS_Mesh cladVectorSelector( cladMeshAdapter );
                auto cladRhsVec = columnRhsVec->select( cladVectorSelector );
                computeTemperatureRhsVector( cladMeshAdapter,
                                             cladTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             cladRhsVec );
            }

            // apply dirichlet rhs correction on f
            if ( bottomPelletBVPOperator ) {
                bottomPelletBVPOperator->modifyRHSvector( columnRhsVec );
            }
            if ( topPelletBVPOperator ) {
                topPelletBVPOperator->modifyRHSvector( columnRhsVec );
            }
            if ( cladBVPOperator ) {
                cladBVPOperator->modifyRHSvector( columnRhsVec );
            }
            {
                auto bottomPelletMat = bottomPelletBVPOperator->getMatrix();
                auto bottomPelletRhs = bottomPelletBVPOperator->subsetOutputVector( columnRhsVec );
                if ( !bottomPelletCor ) {
                    bottomPelletCor = bottomPelletRhs->clone();
                    applyCustomDirichletCondition( bottomPelletRhs,
                                                   bottomPelletCor,
                                                   bottomPelletMeshAdapter,
                                                   bottomPelletConstraints,
                                                   bottomPelletMat );
                } else {
                    applyCustomDirichletCondition( bottomPelletRhs,
                                                   bottomPelletCor,
                                                   bottomPelletMeshAdapter,
                                                   bottomPelletConstraints,
                                                   std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
                } // end if
                AMP_ASSERT( bottomPelletCorptr );

                auto topPelletMat = topPelletBVPOperator->getMatrix();
                auto topPelletRhs = topPelletBVPOperator->subsetOutputVector( columnRhsVec );
                if ( !topPelletCor ) {
                    topPelletCor = topPelletRhs->clone();
                    applyCustomDirichletCondition( topPelletRhs,
                                                   topPelletCor,
                                                   topPelletMeshAdapter,
                                                   topPelletConstraints,
                                                   topPelletMat );
                } else {
                    applyCustomDirichletCondition( topPelletRhs,
                                                   topPelletCor,
                                                   topPelletMeshAdapter,
                                                   topPelletConstraints,
                                                   std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
                } // end if
                AMP_ASSERT( topPelletCorptr );

                auto cladMat = cladBVPOperator->getMatrix();
                auto cladRhs = cladBVPOperator->subsetOutputVector( columnRhsVec );
                if ( !cladCor ) {
                    cladCor = cladRhs->clone();
                    applyCustomDirichletCondition(
                        cladRhs, cladCor, cladMeshAdapter, cladConstraints, cladMat );
                } else {
                    applyCustomDirichletCondition( cladRhs,
                                                   cladCor,
                                                   cladMeshAdapter,
                                                   cladConstraints,
                                                   std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
                } // end if
                AMP_ASSERT( topPelletCorptr );
            }

            // get d
            contactShiftVec->zero();
            bottomPelletTopPelletContactOperator->addShiftToSlave( contactShiftVec );
            bottomPelletCladContactOperator->addShiftToSlave( contactShiftVec );
            topPelletCladContactOperator->addShiftToSlave( contactShiftVec );

            // compute - Kd
            auto rhsCorrectionVec = createVector( dispDofManager, columnVar, split );
            rhsCorrectionVec->zero();
            bottomPelletBVPOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );
            topPelletBVPOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );
            cladBVPOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );

            // f = f - Kd
            columnRhsVec->add( columnRhsVec, rhsCorrectionVec );

            // f^m = f^m + C^T f^s
            // f^s = 0
            bottomPelletTopPelletContactOperator->addSlaveToMaster( columnRhsVec );
            bottomPelletCladContactOperator->addSlaveToMaster( columnRhsVec );
            topPelletCladContactOperator->addSlaveToMaster( columnRhsVec );

            bottomPelletTopPelletContactOperator->setSlaveToZero( columnRhsVec );
            bottomPelletCladContactOperator->setSlaveToZero( columnRhsVec );
            topPelletCladContactOperator->setSlaveToZero( columnRhsVec );

            // u_s = C u_m
            bottomPelletTopPelletContactOperator->copyMasterToSlave( columnSolVec );
            bottomPelletCladContactOperator->copyMasterToSlave( columnSolVec );
            topPelletCladContactOperator->copyMasterToSlave( columnSolVec );

            globalComm.barrier();
            double solveBeginTime = MPI_Wtime();

            linearSolver->apply( columnRhsVec, columnSolVec );

            globalComm.barrier();
            double solveEndTime = MPI_Wtime();
            if ( !rank ) {
                std::cout << "Finished linear solve in " << ( solveEndTime - solveBeginTime )
                          << " seconds." << std::endl;
            }

            // u^s = C u^m + d
            bottomPelletTopPelletContactOperator->copyMasterToSlave( columnSolVec );
            bottomPelletCladContactOperator->copyMasterToSlave( columnSolVec );
            topPelletCladContactOperator->copyMasterToSlave( columnSolVec );

            bottomPelletTopPelletContactOperator->addShiftToSlave( columnSolVec );
            bottomPelletCladContactOperator->addShiftToSlave( columnSolVec );
            topPelletCladContactOperator->addShiftToSlave( columnSolVec );

            computeStressTensor( bottomPelletMeshAdapter,
                                 columnSolVec,
                                 sigma_xx,
                                 sigma_yy,
                                 sigma_zz,
                                 sigma_yz,
                                 sigma_xz,
                                 sigma_xy,
                                 sigma_eff,
                                 fuelMechanicsMaterialModel,
                                 referenceTemperature,
                                 fuelThermalExpansionCoefficient,
                                 tempVec );
            computeStressTensor( topPelletMeshAdapter,
                                 columnSolVec,
                                 sigma_xx,
                                 sigma_yy,
                                 sigma_zz,
                                 sigma_yz,
                                 sigma_xz,
                                 sigma_xy,
                                 sigma_eff,
                                 fuelMechanicsMaterialModel,
                                 referenceTemperature,
                                 fuelThermalExpansionCoefficient,
                                 tempVec );
            computeStressTensor( cladMeshAdapter,
                                 columnSolVec,
                                 sigma_xx,
                                 sigma_yy,
                                 sigma_zz,
                                 sigma_yz,
                                 sigma_xz,
                                 sigma_xy,
                                 sigma_eff,
                                 cladMechanicsMaterialModel,
                                 referenceTemperature,
                                 cladThermalExpansionCoefficient,
                                 tempVec );

            activeSetBeforeUpdateVec->setToScalar( -1.0 );

            const auto &bottomPelletTopPelletActiveSet =
                bottomPelletTopPelletContactOperator->getActiveSet();
            auto bottomPelletTopPelletSizeOfActiveSetBeforeUpdate =
                bottomPelletTopPelletActiveSet.size();

            std::vector<size_t> bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate;
            tempDofManager->getDOFs( bottomPelletTopPelletActiveSet,
                                     bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate );
            AMP_ASSERT( bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate.size() ==
                        bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
            std::vector<double> bottomPelletTopPelletValuesForActiveSet(
                bottomPelletTopPelletSizeOfActiveSetBeforeUpdate, 2.0 );
            activeSetBeforeUpdateVec->setLocalValuesByGlobalID(
                bottomPelletTopPelletSizeOfActiveSetBeforeUpdate,
                &( bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( bottomPelletTopPelletValuesForActiveSet[0] ) );

            std::vector<size_t> bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate;
            dispDofManager->getDOFs( bottomPelletTopPelletActiveSet,
                                     bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate );
            AMP_ASSERT( bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate.size() ==
                        3 * bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
            //
            const auto &bottomPelletCladActiveSet = bottomPelletCladContactOperator->getActiveSet();
            size_t const bottomPelletCladSizeOfActiveSetBeforeUpdate =
                bottomPelletCladActiveSet.size();

            std::vector<size_t> bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate;
            tempDofManager->getDOFs( bottomPelletCladActiveSet,
                                     bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate );
            AMP_ASSERT( bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate.size() ==
                        bottomPelletCladSizeOfActiveSetBeforeUpdate );
            std::vector<double> bottomPelletCladValuesForActiveSet(
                bottomPelletCladSizeOfActiveSetBeforeUpdate, 2.0 );
            activeSetBeforeUpdateVec->setLocalValuesByGlobalID(
                bottomPelletCladSizeOfActiveSetBeforeUpdate,
                &( bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( bottomPelletCladValuesForActiveSet[0] ) );

            std::vector<size_t> bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate;
            dispDofManager->getDOFs( bottomPelletCladActiveSet,
                                     bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate );
            AMP_ASSERT( bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate.size() ==
                        3 * bottomPelletCladSizeOfActiveSetBeforeUpdate );
            //
            const auto &topPelletCladActiveSet = topPelletCladContactOperator->getActiveSet();
            size_t const topPelletCladSizeOfActiveSetBeforeUpdate = topPelletCladActiveSet.size();

            std::vector<size_t> topPelletCladActiveSetTempDOFsIndicesBeforeUpdate;
            tempDofManager->getDOFs( topPelletCladActiveSet,
                                     topPelletCladActiveSetTempDOFsIndicesBeforeUpdate );
            AMP_ASSERT( topPelletCladActiveSetTempDOFsIndicesBeforeUpdate.size() ==
                        topPelletCladSizeOfActiveSetBeforeUpdate );
            std::vector<double> topPelletCladValuesForActiveSet(
                topPelletCladSizeOfActiveSetBeforeUpdate, 2.0 );
            activeSetBeforeUpdateVec->setLocalValuesByGlobalID(
                topPelletCladSizeOfActiveSetBeforeUpdate,
                &( topPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( topPelletCladValuesForActiveSet[0] ) );

            std::vector<size_t> topPelletCladActiveSetDispDOFsIndicesBeforeUpdate;
            dispDofManager->getDOFs( topPelletCladActiveSet,
                                     topPelletCladActiveSetDispDOFsIndicesBeforeUpdate );
            AMP_ASSERT( topPelletCladActiveSetDispDOFsIndicesBeforeUpdate.size() ==
                        3 * topPelletCladSizeOfActiveSetBeforeUpdate );


            {

                if ( scaleSolution != 1.0 ) {
                    columnSolVec->scale( scaleSolution );
                }
                meshAdapter->displaceMesh( columnSolVec );
                columnSolVec->scale( -1.0 );
                meshAdapter->displaceMesh( columnSolVec );
                if ( scaleSolution != 1.0 ) {
                    columnSolVec->scale( -1.0 / scaleSolution );
                } else {
                    columnSolVec->scale( -1.0 );
                }
            }

            size_t nChangesInActiveSet = 0;
            nChangesInActiveSet +=
                bottomPelletTopPelletContactOperator->updateActiveSet( columnSolVec );
            nChangesInActiveSet += bottomPelletCladContactOperator->updateActiveSet( columnSolVec );
            nChangesInActiveSet += topPelletCladContactOperator->updateActiveSet( columnSolVec );

            surfaceTractionVec->zero();
            normalVectorVec->zero();
            contactPressureVec->zero();

            activeSetAfterUpdateVec->setToScalar( -1.0 );

            size_t const bottomPelletTopPelletSizeOfActiveSetAfterUpdate =
                bottomPelletTopPelletActiveSet.size();

            std::vector<size_t> bottomPelletTopPelletActiveSetTempDOFsIndicesAfterUpdate;
            tempDofManager->getDOFs( bottomPelletTopPelletActiveSet,
                                     bottomPelletTopPelletActiveSetTempDOFsIndicesAfterUpdate );
            AMP_ASSERT( bottomPelletTopPelletActiveSetTempDOFsIndicesAfterUpdate.size() ==
                        bottomPelletTopPelletSizeOfActiveSetAfterUpdate );
            std::vector<double> bottomPelletTopPelletValuesForActiveSetAfterUpdate(
                bottomPelletTopPelletSizeOfActiveSetAfterUpdate, 2.0 );
            activeSetAfterUpdateVec->setLocalValuesByGlobalID(
                bottomPelletTopPelletSizeOfActiveSetAfterUpdate,
                &( bottomPelletTopPelletActiveSetTempDOFsIndicesAfterUpdate[0] ),
                &( bottomPelletTopPelletValuesForActiveSetAfterUpdate[0] ) );

            std::vector<size_t> bottomPelletTopPelletActiveSetDispDOFsIndicesAfterUpdate;
            dispDofManager->getDOFs( bottomPelletTopPelletActiveSet,
                                     bottomPelletTopPelletActiveSetDispDOFsIndicesAfterUpdate );
            AMP_ASSERT( bottomPelletTopPelletActiveSetDispDOFsIndicesAfterUpdate.size() ==
                        3 * bottomPelletTopPelletSizeOfActiveSetAfterUpdate );

            std::vector<double> const *bottomPelletTopPelletSlaveVerticesNormalVector;
            std::vector<double> const *bottomPelletTopPelletSlaveVerticesSurfaceTraction;
            bottomPelletTopPelletContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(
                bottomPelletTopPelletSlaveVerticesNormalVector,
                bottomPelletTopPelletSlaveVerticesSurfaceTraction );
            AMP_ASSERT( bottomPelletTopPelletSlaveVerticesSurfaceTraction->size() ==
                        3 * bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
            AMP_ASSERT( bottomPelletTopPelletSlaveVerticesNormalVector->size() ==
                        3 * bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
            surfaceTractionVec->setLocalValuesByGlobalID(
                3 * bottomPelletTopPelletSizeOfActiveSetBeforeUpdate,
                &( bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *bottomPelletTopPelletSlaveVerticesSurfaceTraction )[0] ) );
            normalVectorVec->setLocalValuesByGlobalID(
                3 * bottomPelletTopPelletSizeOfActiveSetBeforeUpdate,
                &( bottomPelletTopPelletActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *bottomPelletTopPelletSlaveVerticesNormalVector )[0] ) );

            std::vector<double> bottomPelletTopPelletSurfaceTractionDOTnormalVector(
                bottomPelletTopPelletSizeOfActiveSetBeforeUpdate );
            for ( size_t kk = 0; kk < bottomPelletTopPelletSizeOfActiveSetBeforeUpdate; ++kk ) {
                bottomPelletTopPelletSurfaceTractionDOTnormalVector[kk] = -compute_scalar_product(
                    &( ( *bottomPelletTopPelletSlaveVerticesSurfaceTraction )[3 * kk] ),
                    &( ( *bottomPelletTopPelletSlaveVerticesNormalVector )[3 * kk] ) );
            } // end for kk
            contactPressureVec->setLocalValuesByGlobalID(
                bottomPelletTopPelletSizeOfActiveSetBeforeUpdate,
                &( bottomPelletTopPelletActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( bottomPelletTopPelletSurfaceTractionDOTnormalVector[0] ) );
            //
            size_t const bottomPelletCladSizeOfActiveSetAfterUpdate =
                bottomPelletCladActiveSet.size();

            std::vector<size_t> bottomPelletCladActiveSetTempDOFsIndicesAfterUpdate;
            tempDofManager->getDOFs( bottomPelletCladActiveSet,
                                     bottomPelletCladActiveSetTempDOFsIndicesAfterUpdate );
            AMP_ASSERT( bottomPelletCladActiveSetTempDOFsIndicesAfterUpdate.size() ==
                        bottomPelletCladSizeOfActiveSetAfterUpdate );
            std::vector<double> bottomPelletCladValuesForActiveSetAfterUpdate(
                bottomPelletCladSizeOfActiveSetAfterUpdate, 2.0 );
            activeSetAfterUpdateVec->setLocalValuesByGlobalID(
                bottomPelletCladSizeOfActiveSetAfterUpdate,
                &( bottomPelletCladActiveSetTempDOFsIndicesAfterUpdate[0] ),
                &( bottomPelletCladValuesForActiveSetAfterUpdate[0] ) );

            std::vector<size_t> bottomPelletCladActiveSetDispDOFsIndicesAfterUpdate;
            dispDofManager->getDOFs( bottomPelletCladActiveSet,
                                     bottomPelletCladActiveSetDispDOFsIndicesAfterUpdate );
            AMP_ASSERT( bottomPelletCladActiveSetDispDOFsIndicesAfterUpdate.size() ==
                        3 * bottomPelletCladSizeOfActiveSetAfterUpdate );

            std::vector<double> const *bottomPelletCladSlaveVerticesNormalVector;
            std::vector<double> const *bottomPelletCladSlaveVerticesSurfaceTraction;
            bottomPelletCladContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(
                bottomPelletCladSlaveVerticesNormalVector,
                bottomPelletCladSlaveVerticesSurfaceTraction );
            AMP_ASSERT( bottomPelletCladSlaveVerticesSurfaceTraction->size() ==
                        3 * bottomPelletCladSizeOfActiveSetBeforeUpdate );
            AMP_ASSERT( bottomPelletCladSlaveVerticesNormalVector->size() ==
                        3 * bottomPelletCladSizeOfActiveSetBeforeUpdate );
            surfaceTractionVec->setLocalValuesByGlobalID(
                3 * bottomPelletCladSizeOfActiveSetBeforeUpdate,
                &( bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *bottomPelletCladSlaveVerticesSurfaceTraction )[0] ) );
            normalVectorVec->setLocalValuesByGlobalID(
                3 * bottomPelletCladSizeOfActiveSetBeforeUpdate,
                &( bottomPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *bottomPelletCladSlaveVerticesNormalVector )[0] ) );

            std::vector<double> bottomPelletCladSurfaceTractionDOTnormalVector(
                bottomPelletCladSizeOfActiveSetBeforeUpdate );
            for ( size_t kk = 0; kk < bottomPelletCladSizeOfActiveSetBeforeUpdate; ++kk ) {
                bottomPelletCladSurfaceTractionDOTnormalVector[kk] = -compute_scalar_product(
                    &( ( *bottomPelletCladSlaveVerticesSurfaceTraction )[3 * kk] ),
                    &( ( *bottomPelletCladSlaveVerticesNormalVector )[3 * kk] ) );
            } // end for kk
            contactPressureVec->setLocalValuesByGlobalID(
                bottomPelletCladSizeOfActiveSetBeforeUpdate,
                &( bottomPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( bottomPelletCladSurfaceTractionDOTnormalVector[0] ) );
            //
            size_t const topPelletCladSizeOfActiveSetAfterUpdate = topPelletCladActiveSet.size();

            std::vector<size_t> topPelletCladActiveSetTempDOFsIndicesAfterUpdate;
            tempDofManager->getDOFs( topPelletCladActiveSet,
                                     topPelletCladActiveSetTempDOFsIndicesAfterUpdate );
            AMP_ASSERT( topPelletCladActiveSetTempDOFsIndicesAfterUpdate.size() ==
                        topPelletCladSizeOfActiveSetAfterUpdate );
            std::vector<double> topPelletCladValuesForActiveSetAfterUpdate(
                topPelletCladSizeOfActiveSetAfterUpdate, 2.0 );
            activeSetAfterUpdateVec->setLocalValuesByGlobalID(
                topPelletCladSizeOfActiveSetAfterUpdate,
                &( topPelletCladActiveSetTempDOFsIndicesAfterUpdate[0] ),
                &( topPelletCladValuesForActiveSetAfterUpdate[0] ) );

            std::vector<size_t> topPelletCladActiveSetDispDOFsIndicesAfterUpdate;
            dispDofManager->getDOFs( topPelletCladActiveSet,
                                     topPelletCladActiveSetDispDOFsIndicesAfterUpdate );
            AMP_ASSERT( topPelletCladActiveSetDispDOFsIndicesAfterUpdate.size() ==
                        3 * topPelletCladSizeOfActiveSetAfterUpdate );

            std::vector<double> const *topPelletCladSlaveVerticesNormalVector;
            std::vector<double> const *topPelletCladSlaveVerticesSurfaceTraction;
            topPelletCladContactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(
                topPelletCladSlaveVerticesNormalVector, topPelletCladSlaveVerticesSurfaceTraction );
            AMP_ASSERT( topPelletCladSlaveVerticesSurfaceTraction->size() ==
                        3 * topPelletCladSizeOfActiveSetBeforeUpdate );
            AMP_ASSERT( topPelletCladSlaveVerticesNormalVector->size() ==
                        3 * topPelletCladSizeOfActiveSetBeforeUpdate );
            surfaceTractionVec->setLocalValuesByGlobalID(
                3 * topPelletCladSizeOfActiveSetBeforeUpdate,
                &( topPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *topPelletCladSlaveVerticesSurfaceTraction )[0] ) );
            normalVectorVec->setLocalValuesByGlobalID(
                3 * topPelletCladSizeOfActiveSetBeforeUpdate,
                &( topPelletCladActiveSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *topPelletCladSlaveVerticesNormalVector )[0] ) );

            std::vector<double> topPelletCladSurfaceTractionDOTnormalVector(
                topPelletCladSizeOfActiveSetBeforeUpdate );
            for ( size_t kk = 0; kk < topPelletCladSizeOfActiveSetBeforeUpdate; ++kk ) {
                topPelletCladSurfaceTractionDOTnormalVector[kk] = -compute_scalar_product(
                    &( ( *topPelletCladSlaveVerticesSurfaceTraction )[3 * kk] ),
                    &( ( *topPelletCladSlaveVerticesNormalVector )[3 * kk] ) );
            } // end for kk
            contactPressureVec->setLocalValuesByGlobalID(
                topPelletCladSizeOfActiveSetBeforeUpdate,
                &( topPelletCladActiveSetTempDOFsIndicesBeforeUpdate[0] ),
                &( topPelletCladSurfaceTractionDOTnormalVector[0] ) );


            if ( scaleSolution != 1.0 ) {
                columnSolVec->scale( scaleSolution );
            }
            meshAdapter->displaceMesh( columnSolVec );
            columnSolVec->scale( -1.0 );
            meshAdapter->displaceMesh( columnSolVec );
            if ( scaleSolution != 1.0 ) {
                columnSolVec->scale( -1.0 / scaleSolution );
            } else {
                columnSolVec->scale( -1.0 );
            }

            if ( !rank ) {
                std::cout << nChangesInActiveSet << " CHANGES IN ACTIVE SET\n";
            }

            if ( nChangesInActiveSet == 0 ) {
                break;
            }
            //    AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
            if ( activeSetIteration == maxActiveSetIterations - 1 ) {
                if ( !rank ) {
                    std::cout << "!!!!!! ACTIVE SET ITERATIONS DID NOT CONVERGE !!!!!!!!\n";
                }
            } // end if
        }     // end for

    } // end for

    meshAdapter->displaceMesh( columnSolVec );

    fout.close();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNodeToGeomType::FaceContactOperator-5" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
