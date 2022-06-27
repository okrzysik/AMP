#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
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

    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->setDecomposition( 1 );

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_driver_" + AMP::Utilities::intToString( rank );
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

    // Load the meshes
    globalComm.barrier();
    double meshBeginTime = MPI_Wtime();

    bool bis                   = input_db->getWithDefault<bool>( "bis", false );
    bool useALittleHelp        = input_db->getWithDefault<bool>( "useALittleHelp", false );
    std::string prefixFileName = input_db->getWithDefault<std::string>( "prefixFileName", "prout" );
    if ( bis ) {
        prefixFileName.append( "bis" );
    }

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( ( bis ? "MeshBis" : "Mesh" ) );
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
    auto columnOperator = std::make_shared<AMP::Operator::ColumnOperator>();

    auto linearSolver_db         = input_db->getDatabase( "LinearSolver" );
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        std::make_shared<AMP::Solver::ColumnSolverParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = columnOperator;
    auto columnPreconditioner =
        std::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    // Get the mechanics material models for the contact operator and for computing stresses
    auto masterModel_db = input_db->getDatabase( "FuelMechanicsMaterialModel" );
    auto masterMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( masterModel_db );
    auto masterMechanicsMaterialModel = std::make_shared<AMP::Operator::IsotropicElasticModel>(
        masterMechanicsMaterialModelParams );
    auto slaveModel_db = input_db->getDatabase(
        ( bis ? "FuelMechanicsMaterialModel" : "CladMechanicsMaterialModel" ) );
    auto slaveMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( slaveModel_db );
    autoslaveMechanicsMaterialModel =
        std::make_shared<AMP::Operator::IsotropicElasticModel>( slaveMechanicsMaterialModelParams );

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    auto contact_db = input_db->getDatabase( ( bis ? "ContactOperatorBis" : "ContactOperator" ) );
    auto contactOperatorParams =
        std::make_shared<AMP::Operator::ContactOperatorParameters>( contact_db );
    contactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    contactOperatorParams->d_DOFManager                   = dispDofManager;
    contactOperatorParams->d_GlobalComm                   = globalComm;
    contactOperatorParams->d_Mesh                         = meshAdapter;
    contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
    contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL
    auto contactOperator = std::make_shared<AMP::Operator::NodeToGeomType::FaceContactOperator>(
        contactOperatorParams );
    contactOperator->initialize();
    contactOperator->setContactIsFrictionless(
        input_db->getWithDefault<bool>( "contactIsFrictionless", false ) );

    bool useML = input_db->getWithDefault<bool>( "useML", false );
    bool cladExpansionConstrained =
        input_db->getWithDefault<bool>( "cladExpansionConstrained", true );
    bool useLevitatingFuel = input_db->getWithDefault<bool>( "useLevitatingFuel", true );
    double scaleSolution   = input_db->getWithDefault<double>( "scaleSolution", 1.0 );
    double shrinkFactor    = input_db->getWithDefault<double>( "shrinkFactor", 0.0 );

    // Build the master and slave operators
    std::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;
    auto masterMeshID      = contactOperator->getMasterMeshID();
    auto masterMeshAdapter = meshAdapter->Subset( masterMeshID );
    //  rotateMesh(masterMeshAdapter);
    if ( masterMeshAdapter.get() != NULL ) {
        auto masterElementPhysicsModel;
        masterBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
        columnOperator->append( masterBVPOperator );

        if ( !useML ) {
            auto masterSolver_db = columnPreconditioner_db->getDatabase( "DummySolver" );
            auto masterSolverParams =
                std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( masterSolver_db );
            masterSolverParams->d_pOperator = masterBVPOperator;
            masterSolverParams->d_comm      = masterMeshAdapter->getComm();
            auto masterSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( masterSolverParams );
            columnPreconditioner->append( masterSolver );
        } else {
            auto masterSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto masterSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( masterSolver_db );
            masterSolverParams->d_pOperator = masterBVPOperator;
            auto masterSolver =
                std::make_shared<AMP::Solver::TrilinosMLSolver>( masterSolverParams );
            columnPreconditioner->append( masterSolver );
        } // end if
    }     // end if

    std::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;
    auto slaveMeshID      = contactOperator->getSlaveMeshID();
    auto slaveMeshAdapter = meshAdapter->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
        slaveBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter,
                ( bis ? "SlaveBVPOperatorBis" : "SlaveBVPOperator" ),
                input_db,
                slaveElementPhysicsModel ) );
        columnOperator->append( slaveBVPOperator );

        if ( !useML ) {
            auto slaveSolver_db = columnPreconditioner_db->getDatabase( "DummySolver" );
            auto slaveSolverParams =
                std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( slaveSolver_db );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            slaveSolverParams->d_comm      = slaveMeshAdapter->getComm();
            auto slaveSolver =
                std::make_shared<AMP::Solver::PetscKrylovSolver>( slaveSolverParams );
            columnPreconditioner->append( slaveSolver );
        } else {
            auto slaveSolver_db = columnPreconditioner_db->getDatabase( "MLSolver" );
            auto slaveSolverParams =
                std::make_shared<AMP::Solver::SolverStrategyParameters>( slaveSolver_db );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
            auto slaveSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( slaveSolverParams );
            columnPreconditioner->append( slaveSolver );
        } // end if
    }     // end if

    auto contactPreconditioner_db = columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    auto contactPreconditionerParams =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
            contactPreconditioner_db );
    contactPreconditionerParams->d_pOperator = contactOperator;
    autocontactPreconditioner =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolver>( contactPreconditionerParams );
    columnPreconditioner->append( contactPreconditioner );

    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> masterConstraints;
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> slaveConstraints;
    if ( bis ) {
        if ( !useLevitatingFuel ) {
            double fuelOuterRadius = input_db->getScalar<double>( "FuelOuterRadius" );
            makeConstraintsOnFuel(
                masterMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
                fuelOuterRadius,
                masterConstraints,
                true );
            makeConstraintsOnFuel(
                masterMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
                fuelOuterRadius,
                masterConstraints,
                false );
            makeConstraintsOnFuel(
                slaveMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
                fuelOuterRadius,
                slaveConstraints,
                true );
            makeConstraintsOnFuel(
                slaveMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
                fuelOuterRadius,
                slaveConstraints,
                false );
        }
    } else {
        if ( !cladExpansionConstrained ) {
            double cladInnerRadius = input_db->getScalar<double>( "CladInnerRadius" );
            double cladOuterRadius = input_db->getScalar<double>( "CladOuterRadius" );
            double cladHeight      = input_db->getScalar<double>( "CladHeight" );
            makeConstraintsOnClad( slaveMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex ),
                                   cladInnerRadius,
                                   cladOuterRadius,
                                   cladHeight,
                                   slaveConstraints );
        }
        if ( !useLevitatingFuel ) {
            double fuelOuterRadius = input_db->getScalar<double>( "FuelOuterRadius" );
            makeConstraintsOnFuel(
                masterMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
                fuelOuterRadius,
                masterConstraints,
                false );
            makeConstraintsOnFuel(
                masterMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
                fuelOuterRadius,
                masterConstraints,
                true );
        }
    } // end if

    // Items for computing the RHS correction due to thermal expansion
    auto masterTemperatureRhs_db = input_db->getDatabase( "MasterTemperatureRHSVectorCorrection" );
    auto slaveTemperatureRhs_db  = input_db->getDatabase(
        ( bis ? "MasterTemperatureRHSVectorCorrection" : "SlaveTemperatureRHSVectorCorrection" ) );
    auto tempVar        = std::make_shared<AMP::LinearAlgebra::Variable>( "temperature" );
    auto dispVar        = columnOperator->getOutputVariable();
    auto tempDofManager = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, 1, split );
    auto tempVec    = AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    auto refTempVec = tempVec->cloneVector();

    auto sigma_xx = AMP::LinearAlgebra::createVector(
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

    computeFuelTemperature( masterMeshAdapter,
                            tempVec,
                            fuelOuterRadius,
                            fuelOuterRadiusTemperature,
                            linearHeatGenerationRate,
                            dummyFuelThermalConductivity );
    if ( !bis ) {
        computeCladTemperature( slaveMeshAdapter,
                                tempVec,
                                cladInnerRadius,
                                cladOuterRadius,
                                linearHeatGenerationRate,
                                cladOuterRadiusTemperature,
                                cladThermalConductivity );
    } else {
        computeFuelTemperature( slaveMeshAdapter,
                                tempVec,
                                fuelOuterRadius,
                                fuelOuterRadiusTemperature,
                                linearHeatGenerationRate,
                                dummyFuelThermalConductivity );
    } // end if

    globalComm.barrier();
    double tempCompEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished computing the temperature profile in "
                  << ( tempCompEndTime - tempCompBeginTime ) << " seconds." << std::endl;
    }

    // AMP::LinearAlgebra::VS_Mesh slaveVectorSelector(slaveMeshAdapter);
    // auto slaveTempVec = tempVec->select(slaveVectorSelector, tempVar->getName());
    // slaveTempVec->setToScalar(900.0);
    // auto tmp_db = masterTemperatureRhs_db->getDatabase("RhsMaterialModel");
    // double masterThermalExpansionCoefficient =
    //    tmp_db->getScalar<double>("THERMAL_EXPANSION_COEFFICIENT");
    double masterThermalExpansionCoefficient =
        ( masterTemperatureRhs_db->getDatabase( "RhsMaterialModel" ) )
            ->getScalar<double>( "THERMAL_EXPANSION_COEFFICIENT" );
    double slaveThermalExpansionCoefficient =
        ( slaveTemperatureRhs_db->getDatabase( "RhsMaterialModel" ) )
            ->getScalar<double>( "THERMAL_EXPANSION_COEFFICIENT" );
    contactOperator->uglyHack(
        tempVec, tempDofManager, masterThermalExpansionCoefficient, referenceTemperature );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnVar    = columnOperator->getOutputVariable();
    auto columnSolVec = AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    auto columnRhsVec = AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();
    AMP::LinearAlgebra::Vector::shared_ptr masterCor;
    AMP::LinearAlgebra::Vector::shared_ptr slaveCor;
    auto activeSetBeforeUpdateVec = sigma_eff->cloneVector();
    auto activeSetAfterUpdateVec  = sigma_eff->cloneVector();
    auto contactPressureVec       = sigma_eff->cloneVector();
    auto surfaceTractionVec       = columnSolVec->cloneVector();
    auto normalVectorVec          = columnSolVec->cloneVector();

    if ( ( !bis ) && ( shrinkFactor != 0.0 ) ) {
        AMP_ASSERT( ( shrinkFactor > 0.0 ) && ( shrinkFactor < 1.0 ) );
        shrinkMesh( slaveMeshAdapter, shrinkFactor );
    }

    if ( bis && useALittleHelp ) {
        // auto zDispVec = columnSolVec->select(AMP::LinearAlgebra::VS_Stride(2,3), "help");
        auto it       = slaveMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 );
        auto it_begin = it.begin();
        auto it_end   = it.end();
        std::vector<double> coord;
        std::vector<size_t> dofs;
        double radius = -1.0;
        for ( it = it_begin; it != it_end; ++it ) {
            coord  = it->coord();
            radius = std::sqrt( std::pow( coord[0], 2 ) + std::pow( coord[1], 2 ) );
            if ( radius > 0.0022 ) {
                std::cout << radius << "  " << coord[2] << "\n";
                dispDofManager->getDOFs( it->globalID(), dofs );
                //        zDispVec->setValueByGlobalID(dofs[2], 0.00005);
                columnSolVec->setValueByGlobalID( dofs[2], 0.00005 );
            } // end if
        }     // end for
        contactOperator->updateActiveSetWithALittleHelp( columnSolVec );
        //    zDispVec->zero();
        columnSolVec->zero();
    } else {
        bool skipDisplaceMesh = true;
        contactOperator->updateActiveSet( nullVec, skipDisplaceMesh );
    } // end if

    auto contactShiftVec = createVector( dispDofManager, columnVar, split );
    contactShiftVec->zero();

    siloWriter->registerVector(
        columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    siloWriter->registerVector( tempVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Temperature" );
    siloWriter->registerVector( sigma_eff, meshAdapter, AMP::Mesh::GeomType::Vertex, "vonMises" );
    siloWriter->registerVector( sigma_xx, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xx" );
    siloWriter->registerVector( sigma_yy, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_yy" );
    siloWriter->registerVector( sigma_zz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_zz" );
    siloWriter->registerVector( sigma_yz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_yz" );
    siloWriter->registerVector( sigma_xz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xz" );
    siloWriter->registerVector( sigma_xy, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xy" );
    siloWriter->registerVector( activeSetBeforeUpdateVec,
                                meshAdapter,
                                AMP::Mesh::GeomType::Vertex,
                                "ActiveSetBeforeUpdate" );
    siloWriter->registerVector(
        activeSetAfterUpdateVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ActiveSetAfterUpdate" );
    siloWriter->registerVector(
        surfaceTractionVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Traction" );
    siloWriter->registerVector(
        normalVectorVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Normal" );
    siloWriter->registerVector(
        contactPressureVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ContactPressure" );
    siloWriter->registerVector(
        contactShiftVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Shift" );
    siloWriter->writeFile( prefixFileName.c_str(), 0 );

    columnSolVec->zero();
    columnOperator->append( contactOperator );

    // Build a matrix shell operator to use the column operator with the petsc krylov solvers
    auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
    auto matrixShellParams =
        std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase );
    auto matrixShellOperator =
        std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );

    int numMasterLocalNodes = 0;
    int numSlaveLocalNodes  = 0;
    if ( masterMeshAdapter.get() != NULL ) {
        numMasterLocalNodes = masterMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( slaveMeshAdapter.get() != NULL ) {
        numSlaveLocalNodes = slaveMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    int matLocalSize = dofsPerNode * ( numMasterLocalNodes + numSlaveLocalNodes );
    AMP_ASSERT( matLocalSize == static_cast<int>( dispDofManager->numLocalDOF() ) );
    matrixShellOperator->setComm( globalComm );
    matrixShellOperator->setMatLocalRowSize( matLocalSize );
    matrixShellOperator->setMatLocalColumnSize( matLocalSize );
    matrixShellOperator->setOperator( columnOperator );

    auto linearSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = matrixShellOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = columnPreconditioner;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    auto fullThermalLoadingTempMinusRefTempVec = tempVec->cloneVector();
    fullThermalLoadingTempMinusRefTempVec->subtract( tempVec, refTempVec );

    int TOTO_count = 0;
    size_t const maxThermalLoadingIterations =
        input_db->getWithDefault<size_t>( "maxThermalLoadingIterations", 5 );
    for ( size_t thermalLoadingIteration = 0; thermalLoadingIteration < maxThermalLoadingIterations;
          ++thermalLoadingIteration ) {
        if ( !rank ) {
            std::cout << "THERMAL LOADING " << thermalLoadingIteration + 1 << "/"
                      << maxThermalLoadingIterations << "\n";
        }
        double scalingFactor = static_cast<double>( thermalLoadingIteration + 1 ) /
                               static_cast<double>( maxThermalLoadingIterations );
        tempVec->axpy( scalingFactor, fullThermalLoadingTempMinusRefTempVec, refTempVec );

        size_t const maxActiveSetIterations =
            input_db->getWithDefault<size_t>( "maxActiveSetIterations", 5 );
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
                AMP::LinearAlgebra::VS_Mesh slaveVectorSelector( slaveMeshAdapter );
                auto slaveRhsVec = columnRhsVec->select( slaveVectorSelector, dispVar->getName() );
                computeTemperatureRhsVector( slaveMeshAdapter,
                                             slaveTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             slaveRhsVec );
            }
            {
                AMP::LinearAlgebra::VS_Mesh masterVectorSelector( masterMeshAdapter );
                auto masterRhsVec =
                    columnRhsVec->select( masterVectorSelector, dispVar->getName() );
                computeTemperatureRhsVector( masterMeshAdapter,
                                             masterTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             masterRhsVec );
            }

            // apply dirichlet rhs correction on f
            if ( masterBVPOperator.get() != NULL ) {
                masterBVPOperator->modifyRHSvector( columnRhsVec );
            } // end if
            if ( slaveBVPOperator.get() != NULL ) {
                slaveBVPOperator->modifyRHSvector( columnRhsVec );
            } // end if

            {
                auto masterMat = masterBVPOperator->getMatrix();
                auto masterRhs = masterBVPOperator->subsetOutputVector( columnRhsVec );
                if ( masterCor.get() == nullptr ) {
                    masterCor = masterRhs->cloneVector();
                    applyCustomDirichletCondition(
                        masterRhs, masterCor, meshAdapter, masterConstraints, masterMat );
                } else {
                    applyCustomDirichletCondition( masterRhs,
                                                   masterCor,
                                                   meshAdapter,
                                                   masterConstraints,
                                                   std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
                } // end if
                AMP_ASSERT( masterCor.get() != nullptr );
            }
            {
                auto slaveMat = slaveBVPOperator->getMatrix();
                auto slaveRhs = slaveBVPOperator->subsetOutputVector( columnRhsVec );
                if ( slaveCor.get() == nullptr ) {
                    slaveCor = slaveRhs->cloneVector();
                    applyCustomDirichletCondition(
                        slaveRhs, slaveCor, meshAdapter, slaveConstraints, slaveMat );
                } else {
                    applyCustomDirichletCondition( slaveRhs,
                                                   slaveCor,
                                                   meshAdapter,
                                                   slaveConstraints,
                                                   std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
                } // end if
                AMP_ASSERT( slaveCor.get() != nullptr );
            }

            // get d
            // auto contactShiftVec = createVector(dispDofManager, columnVar, split);
            contactShiftVec->zero();
            contactOperator->addShiftToSlave( contactShiftVec );
            //  contactOperator->addShiftToSlave(columnSolVec);

            // compute - Kd
            auto rhsCorrectionVec = createVector( dispDofManager, columnVar, split );
            rhsCorrectionVec->zero();
            masterBVPOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );
            slaveBVPOperator->apply( nullVec, contactShiftVec, rhsCorrectionVec, -1.0, 0.0 );
            //  columnOperator->apply(nullVec, columnSolVec, rhsCorrectionVec, -1.0, 0.0);
            //  columnOperator->append(contactOperator);

            // f = f - Kd
            columnRhsVec->add( columnRhsVec, rhsCorrectionVec );

            // f^m = f^m + C^T f^s
            // f^s = 0
            contactOperator->addSlaveToMaster( columnRhsVec );
            contactOperator->setSlaveToZero( columnRhsVec );

            // u_s = C u_m
            contactOperator->copyMasterToSlave( columnSolVec );

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
            contactOperator->copyMasterToSlave( columnSolVec );
            contactOperator->addShiftToSlave( columnSolVec );

            computeStressTensor( masterMeshAdapter,
                                 columnSolVec,
                                 sigma_xx,
                                 sigma_yy,
                                 sigma_zz,
                                 sigma_yz,
                                 sigma_xz,
                                 sigma_xy,
                                 sigma_eff,
                                 masterMechanicsMaterialModel,
                                 referenceTemperature,
                                 masterThermalExpansionCoefficient,
                                 tempVec );
            computeStressTensor( slaveMeshAdapter,
                                 columnSolVec,
                                 sigma_xx,
                                 sigma_yy,
                                 sigma_zz,
                                 sigma_yz,
                                 sigma_xz,
                                 sigma_xy,
                                 sigma_eff,
                                 slaveMechanicsMaterialModel,
                                 referenceTemperature,
                                 slaveThermalExpansionCoefficient,
                                 tempVec );

            const auto &activeSet                    = contactOperator->getActiveSet();
            size_t const sizeOfActiveSetBeforeUpdate = activeSet.size();

            std::vector<size_t> activeSetTempDOFsIndicesBeforeUpdate;
            tempDofManager->getDOFs( activeSet, activeSetTempDOFsIndicesBeforeUpdate );
            AMP_ASSERT( activeSetTempDOFsIndicesBeforeUpdate.size() ==
                        sizeOfActiveSetBeforeUpdate );
            std::vector<double> valuesForActiveSetBeforeUpdate( sizeOfActiveSetBeforeUpdate, 2.0 );
            activeSetBeforeUpdateVec->setToScalar( -1.0 );
            activeSetBeforeUpdateVec->setLocalValuesByGlobalID(
                sizeOfActiveSetBeforeUpdate,
                &( activeSetTempDOFsIndicesBeforeUpdate[0] ),
                &( valuesForActiveSetBeforeUpdate[0] ) );

            std::vector<size_t> activeSetDispDOFsIndicesBeforeUpdate;
            dispDofManager->getDOFs( activeSet, activeSetDispDOFsIndicesBeforeUpdate );
            AMP_ASSERT( activeSetDispDOFsIndicesBeforeUpdate.size() ==
                        3 * sizeOfActiveSetBeforeUpdate );


            if ( scaleSolution != 1.0 ) {
                columnSolVec->scale( scaleSolution );
            }
            meshAdapter->displaceMesh( columnSolVec );
            siloWriter->writeFile( prefixFileName.c_str(), TOTO_count );
            columnSolVec->scale( -1.0 );
            meshAdapter->displaceMesh( columnSolVec );
            if ( scaleSolution != 1.0 ) {
                columnSolVec->scale( -1.0 / scaleSolution );
            } else {
                columnSolVec->scale( -1.0 );
            }

            size_t nChangesInActiveSet = contactOperator->updateActiveSet( columnSolVec );

            size_t const sizeOfActiveSetAfterUpdate = activeSet.size();

            std::vector<size_t> activeSetTempDOFsIndicesAfterUpdate;
            tempDofManager->getDOFs( activeSet, activeSetTempDOFsIndicesAfterUpdate );
            AMP_ASSERT( activeSetTempDOFsIndicesAfterUpdate.size() == sizeOfActiveSetAfterUpdate );
            std::vector<double> valuesForActiveSetAfterUpdate( sizeOfActiveSetAfterUpdate, 2.0 );
            activeSetAfterUpdateVec->setToScalar( -1.0 );
            activeSetAfterUpdateVec->setLocalValuesByGlobalID(
                sizeOfActiveSetAfterUpdate,
                &( activeSetTempDOFsIndicesAfterUpdate[0] ),
                &( valuesForActiveSetAfterUpdate[0] ) );

            std::vector<size_t> activeSetDispDOFsIndicesAfterUpdate;
            dispDofManager->getDOFs( activeSet, activeSetDispDOFsIndicesAfterUpdate );
            AMP_ASSERT( activeSetDispDOFsIndicesAfterUpdate.size() ==
                        3 * sizeOfActiveSetAfterUpdate );

            std::vector<double> const *slaveVerticesNormalVectorBeforeUpdate;
            std::vector<double> const *slaveVerticesSurfaceTractionBeforeUpdate;
            contactOperator->getSlaveVerticesNormalVectorAndSurfaceTraction(
                slaveVerticesNormalVectorBeforeUpdate, slaveVerticesSurfaceTractionBeforeUpdate );
            AMP_ASSERT( slaveVerticesSurfaceTractionBeforeUpdate->size() ==
                        3 * sizeOfActiveSetBeforeUpdate );
            AMP_ASSERT( slaveVerticesNormalVectorBeforeUpdate->size() ==
                        3 * sizeOfActiveSetBeforeUpdate );
            surfaceTractionVec->zero();
            surfaceTractionVec->setLocalValuesByGlobalID(
                3 * sizeOfActiveSetBeforeUpdate,
                &( activeSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *slaveVerticesSurfaceTractionBeforeUpdate )[0] ) );
            normalVectorVec->zero();
            normalVectorVec->setLocalValuesByGlobalID(
                3 * sizeOfActiveSetBeforeUpdate,
                &( activeSetDispDOFsIndicesBeforeUpdate[0] ),
                &( ( *slaveVerticesNormalVectorBeforeUpdate )[0] ) );

            std::vector<double> surfaceTractionDOTnormalVector( sizeOfActiveSetBeforeUpdate );
            for ( size_t kk = 0; kk < sizeOfActiveSetBeforeUpdate; ++kk ) {
                surfaceTractionDOTnormalVector[kk] = -compute_scalar_product(
                    &( ( *slaveVerticesSurfaceTractionBeforeUpdate )[3 * kk] ),
                    &( ( *slaveVerticesNormalVectorBeforeUpdate )[3 * kk] ) );
            } // end for kk
            contactPressureVec->zero();
            contactPressureVec->setLocalValuesByGlobalID(
                sizeOfActiveSetBeforeUpdate,
                &( activeSetTempDOFsIndicesBeforeUpdate[0] ),
                &( surfaceTractionDOTnormalVector[0] ) );

            if ( scaleSolution != 1.0 ) {
                columnSolVec->scale( scaleSolution );
            }
            meshAdapter->displaceMesh( columnSolVec );
            siloWriter->writeFile( prefixFileName.c_str(), TOTO_count );
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
    exeNames.push_back( "testNodeToGeomType::FaceContactOperator-4" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
