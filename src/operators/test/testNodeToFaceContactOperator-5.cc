
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "externVars.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/Writer.h"

#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/TrilinosMatrixShellOperator.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/contact/NodeToGeomType::FaceContactOperator.h"
#include "operators/mechanics/ConstructLinearMechanicsRHSVector.h"
#include "operators/mechanics/IsotropicElasticModel.h"
#include "operators/mechanics/MechanicsLinearFEOperator.h"
#include "operators/mechanics/MechanicsMaterialModel.h"
#include "operators/mechanics/MechanicsModelParameters.h"
#include "operators/petsc/PetscMatrixShellOperator.h"

#include "solvers/ColumnSolver.h"
#include "solvers/ConstraintsEliminationSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

#include "ampmesh/euclidean_geometry_tools.h"
#include "ampmesh/latex_visualization_tools.h"
#include <fstream>
#include <set>

#include "testNodeToGeomType::FaceContactOperator.h"


void myTest( AMP::UnitTest *ut, std::string exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->setDecomposition( 1 );
#endif

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    std::fstream fout;
    std::string fileName = "debug_driver_" + AMP::Utilities::intToString( rank );
    fout.open( fileName.c_str(), std::fstream::out );

    // Load the input file
    globalComm.barrier();
    double inpReadBeginTime = MPI_Wtime();

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    globalComm.barrier();
    double inpReadEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished parsing the input file in " << ( inpReadEndTime - inpReadBeginTime )
                  << " seconds." << std::endl;
    }

    bool useML = input_db->getBoolWithDefault( "useML", false );
    bool cladExpansionConstrained =
        input_db->getBoolWithDefault( "cladExpansionConstrained", true );
    bool useLevitatingFuel     = input_db->getBoolWithDefault( "useLevitatingFuel", true );
    std::string prefixFileName = input_db->getStringWithDefault( "prefixFileName", "TATA_0" );
    double scaleSolution       = input_db->getDoubleWithDefault( "scaleSolution", 1.0 );
    double cladNeedALittleHelp = input_db->getDoubleWithDefault( "cladNeedALittleHelp", 0.0 );
    double fuelNeedALittleHelp = input_db->getDoubleWithDefault( "fuelNeedALittleHelp", -1.0 );
    bool contactIsFrictionless = input_db->getBoolWithDefault( "contactIsFrictionless", false );
    double shrinkFactor        = input_db->getDoubleWithDefault( "shrinkFactor", 0.0 );

    // Load the meshes
    globalComm.barrier();
    double meshBeginTime = MPI_Wtime();

    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    meshParams->setComm( globalComm );
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

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
    AMP::Discretization::DOFManager::shared_ptr dispDofManager =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, dofsPerNode, split );

    // Build a column operator and a column preconditioner
    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(
        new AMP::Operator::ColumnOperator( emptyParams ) );

    AMP::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );
    AMP::shared_ptr<AMP::Database> columnPreconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::ColumnSolverParameters> columnPreconditionerParams(
        new AMP::Solver::ColumnSolverParameters( columnPreconditioner_db ) );
    columnPreconditionerParams->d_pOperator = columnOperator;
    AMP::shared_ptr<AMP::Solver::ColumnSolver> columnPreconditioner(
        new AMP::Solver::ColumnSolver( columnPreconditionerParams ) );

    // Get the mechanics material models for the contact operator and for computing stresses
    AMP::shared_ptr<AMP::Database> fuelModel_db =
        input_db->getDatabase( "FuelMechanicsMaterialModel" );
    AMP::shared_ptr<AMP::Operator::MechanicsModelParameters> fuelMechanicsMaterialModelParams(
        new AMP::Operator::MechanicsModelParameters( fuelModel_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> fuelMechanicsMaterialModel(
        new AMP::Operator::IsotropicElasticModel( fuelMechanicsMaterialModelParams ) );
    AMP::shared_ptr<AMP::Database> cladModel_db =
        input_db->getDatabase( "CladMechanicsMaterialModel" );
    AMP::shared_ptr<AMP::Operator::MechanicsModelParameters> cladMechanicsMaterialModelParams(
        new AMP::Operator::MechanicsModelParameters( cladModel_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> cladMechanicsMaterialModel(
        new AMP::Operator::IsotropicElasticModel( cladMechanicsMaterialModelParams ) );

    // Build the contact operators
    AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator> bottomPelletTopPelletContactOperator;
    AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator> bottomPelletCladContactOperator;
    AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator> topPelletCladContactOperator;

    AMP::shared_ptr<AMP::Database> bottomPelletTopPelletContact_db =
        input_db->getDatabase( "BottomPelletTopPelletContactOperator" );
    AMP::shared_ptr<AMP::Operator::ContactOperatorParameters>
        bottomPelletTopPelletContactOperatorParams(
            new AMP::Operator::ContactOperatorParameters( bottomPelletTopPelletContact_db ) );
    bottomPelletTopPelletContactOperatorParams->d_DOFsPerNode = dofsPerNode;
    bottomPelletTopPelletContactOperatorParams->d_DOFManager  = dispDofManager;
    bottomPelletTopPelletContactOperatorParams->d_GlobalComm  = globalComm;
    bottomPelletTopPelletContactOperatorParams->d_Mesh        = meshAdapter;
    bottomPelletTopPelletContactOperatorParams->d_MasterMechanicsMaterialModel =
        fuelMechanicsMaterialModel;
    bottomPelletTopPelletContactOperatorParams->reset();
    bottomPelletTopPelletContactOperator =
        AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator>(
            new AMP::Operator::NodeToGeomType::FaceContactOperator(
                bottomPelletTopPelletContactOperatorParams ) );
    bottomPelletTopPelletContactOperator->initialize();
    bottomPelletTopPelletContactOperator->setContactIsFrictionless( contactIsFrictionless );
    //  bottomPelletTopPelletContactOperator->setContactIsFrictionless(bottomPelletTopPelletContact_db->getBoolWithDefault("ContactIsFrictionless",
    //  false));

    AMP::shared_ptr<AMP::Database> bottomPelletCladContact_db =
        input_db->getDatabase( "BottomPelletCladContactOperator" );
    AMP::shared_ptr<AMP::Operator::ContactOperatorParameters> bottomPelletCladContactOperatorParams(
        new AMP::Operator::ContactOperatorParameters( bottomPelletCladContact_db ) );
    bottomPelletCladContactOperatorParams->d_DOFsPerNode = dofsPerNode;
    bottomPelletCladContactOperatorParams->d_DOFManager  = dispDofManager;
    bottomPelletCladContactOperatorParams->d_GlobalComm  = globalComm;
    bottomPelletCladContactOperatorParams->d_Mesh        = meshAdapter;
    bottomPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel =
        fuelMechanicsMaterialModel;
    bottomPelletCladContactOperatorParams->reset();
    bottomPelletCladContactOperator = AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator>(
        new AMP::Operator::NodeToGeomType::FaceContactOperator( bottomPelletCladContactOperatorParams ) );
    bottomPelletCladContactOperator->initialize();
    bottomPelletCladContactOperator->setContactIsFrictionless( contactIsFrictionless );
    //  bottomPelletCladContactOperator->setContactIsFrictionless(bottomPelletCladContact_db->getBoolWithDefault("ContactIsFrictionless",
    //  false));

    AMP::shared_ptr<AMP::Database> topPelletCladContact_db =
        input_db->getDatabase( "TopPelletCladContactOperator" );
    AMP::shared_ptr<AMP::Operator::ContactOperatorParameters> topPelletCladContactOperatorParams(
        new AMP::Operator::ContactOperatorParameters( topPelletCladContact_db ) );
    topPelletCladContactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    topPelletCladContactOperatorParams->d_DOFManager                   = dispDofManager;
    topPelletCladContactOperatorParams->d_GlobalComm                   = globalComm;
    topPelletCladContactOperatorParams->d_Mesh                         = meshAdapter;
    topPelletCladContactOperatorParams->d_MasterMechanicsMaterialModel = fuelMechanicsMaterialModel;
    topPelletCladContactOperatorParams->reset();
    topPelletCladContactOperator = AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator>(
        new AMP::Operator::NodeToGeomType::FaceContactOperator( topPelletCladContactOperatorParams ) );
    topPelletCladContactOperator->initialize();
    topPelletCladContactOperator->setContactIsFrictionless( contactIsFrictionless );
    //  topPelletCladContactOperator->setContactIsFrictionless(topPelletCladContact_db->getBoolWithDefault("ContactIsFrictionless",
    //  false));

    // Build the BVP operators
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bottomPelletBVPOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> topPelletBVPOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> cladBVPOperator;


    AMP::Mesh::MeshID bottomPelletMeshID = bottomPelletTopPelletContactOperator->getMasterMeshID();
    AMP_ASSERT( bottomPelletMeshID == bottomPelletCladContactOperator->getMasterMeshID() );
    AMP::Mesh::Mesh::shared_ptr bottomPelletMeshAdapter = meshAdapter->Subset( bottomPelletMeshID );
    if ( bottomPelletMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> bottomPelletElementPhysicsModel;
        bottomPelletBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( bottomPelletMeshAdapter,
                                                            "BottomPelletBVPOperator",
                                                            input_db,
                                                            bottomPelletElementPhysicsModel ) );
        columnOperator->append( bottomPelletBVPOperator );

        if ( !useML ) {
            AMP::shared_ptr<AMP::Database> bottomPelletSolver_db =
                columnPreconditioner_db->getDatabase( "DummySolver" );
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> bottomPelletSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( bottomPelletSolver_db ) );
            bottomPelletSolverParams->d_pOperator = bottomPelletBVPOperator;
            bottomPelletSolverParams->d_comm      = bottomPelletMeshAdapter->getComm();
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> bottomPelletSolver(
                new AMP::Solver::PetscKrylovSolver( bottomPelletSolverParams ) );
            columnPreconditioner->append( bottomPelletSolver );
        } else {
            AMP::shared_ptr<AMP::Database> bottomPelletSolver_db =
                columnPreconditioner_db->getDatabase( "MLSolver" );
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> bottomPelletSolverParams(
                new AMP::Solver::SolverStrategyParameters( bottomPelletSolver_db ) );
            bottomPelletSolverParams->d_pOperator = bottomPelletBVPOperator;
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> bottomPelletSolver(
                new AMP::Solver::TrilinosMLSolver( bottomPelletSolverParams ) );
            columnPreconditioner->append( bottomPelletSolver );
        } // end if
    }     // end if

    AMP::Mesh::MeshID topPelletMeshID = bottomPelletTopPelletContactOperator->getSlaveMeshID();
    AMP_ASSERT( topPelletMeshID == topPelletCladContactOperator->getMasterMeshID() );
    AMP::Mesh::Mesh::shared_ptr topPelletMeshAdapter = meshAdapter->Subset( topPelletMeshID );
    if ( topPelletMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> topPelletElementPhysicsModel;
        topPelletBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator( topPelletMeshAdapter,
                                                            "TopPelletBVPOperator",
                                                            input_db,
                                                            topPelletElementPhysicsModel ) );
        columnOperator->append( topPelletBVPOperator );

        if ( !useML ) {
            AMP::shared_ptr<AMP::Database> topPelletSolver_db =
                columnPreconditioner_db->getDatabase( "DummySolver" );
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> topPelletSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( topPelletSolver_db ) );
            topPelletSolverParams->d_pOperator = topPelletBVPOperator;
            topPelletSolverParams->d_comm      = topPelletMeshAdapter->getComm();
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> topPelletSolver(
                new AMP::Solver::PetscKrylovSolver( topPelletSolverParams ) );
            columnPreconditioner->append( topPelletSolver );
        } else {
            AMP::shared_ptr<AMP::Database> topPelletSolver_db =
                columnPreconditioner_db->getDatabase( "MLSolver" );
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> topPelletSolverParams(
                new AMP::Solver::SolverStrategyParameters( topPelletSolver_db ) );
            topPelletSolverParams->d_pOperator = topPelletBVPOperator;
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> topPelletSolver(
                new AMP::Solver::TrilinosMLSolver( topPelletSolverParams ) );
            columnPreconditioner->append( topPelletSolver );
        } // end if
    }     // end if

    AMP::Mesh::MeshID cladMeshID = bottomPelletCladContactOperator->getSlaveMeshID();
    AMP_ASSERT( cladMeshID == topPelletCladContactOperator->getSlaveMeshID() );
    AMP::Mesh::Mesh::shared_ptr cladMeshAdapter = meshAdapter->Subset( cladMeshID );
    if ( cladMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> cladElementPhysicsModel;
        cladBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                cladMeshAdapter, "CladBVPOperator", input_db, cladElementPhysicsModel ) );
        columnOperator->append( cladBVPOperator );

        if ( useML ) {
            AMP::shared_ptr<AMP::Database> cladSolver_db =
                columnPreconditioner_db->getDatabase( "MLSolver" );
            AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> cladSolverParams(
                new AMP::Solver::SolverStrategyParameters( cladSolver_db ) );
            cladSolverParams->d_pOperator = cladBVPOperator;
            AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> cladSolver(
                new AMP::Solver::TrilinosMLSolver( cladSolverParams ) );
            columnPreconditioner->append( cladSolver );
        } else {
            AMP::shared_ptr<AMP::Database> cladSolver_db =
                columnPreconditioner_db->getDatabase( "DummySolver" );
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> cladSolverParams(
                new AMP::Solver::PetscKrylovSolverParameters( cladSolver_db ) );
            cladSolverParams->d_pOperator = cladBVPOperator;
            cladSolverParams->d_comm      = cladMeshAdapter->getComm();
            AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> cladSolver(
                new AMP::Solver::PetscKrylovSolver( cladSolverParams ) );
            columnPreconditioner->append( cladSolver );
        } // end if

    } // end if


    {
        AMP::shared_ptr<AMP::Database> bottomPelletTopPelletContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters>
            bottomPelletTopPelletContactPreconditionerParams(
                new AMP::Solver::ConstraintsEliminationSolverParameters(
                    bottomPelletTopPelletContactPreconditioner_db ) );
        bottomPelletTopPelletContactPreconditionerParams->d_pOperator =
            bottomPelletTopPelletContactOperator;
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver>
            bottomPelletTopPelletContactPreconditioner(
                new AMP::Solver::ConstraintsEliminationSolver(
                    bottomPelletTopPelletContactPreconditionerParams ) );
        columnPreconditioner->append( bottomPelletTopPelletContactPreconditioner );
    }
    {
        AMP::shared_ptr<AMP::Database> bottomPelletCladContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters>
            bottomPelletCladContactPreconditionerParams(
                new AMP::Solver::ConstraintsEliminationSolverParameters(
                    bottomPelletCladContactPreconditioner_db ) );
        bottomPelletCladContactPreconditionerParams->d_pOperator = bottomPelletCladContactOperator;
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver>
            bottomPelletCladContactPreconditioner( new AMP::Solver::ConstraintsEliminationSolver(
                bottomPelletCladContactPreconditionerParams ) );
        columnPreconditioner->append( bottomPelletCladContactPreconditioner );
    }
    {
        AMP::shared_ptr<AMP::Database> topPelletCladContactPreconditioner_db =
            columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters>
            topPelletCladContactPreconditionerParams(
                new AMP::Solver::ConstraintsEliminationSolverParameters(
                    topPelletCladContactPreconditioner_db ) );
        topPelletCladContactPreconditionerParams->d_pOperator = topPelletCladContactOperator;
        AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver>
            topPelletCladContactPreconditioner( new AMP::Solver::ConstraintsEliminationSolver(
                topPelletCladContactPreconditionerParams ) );
        columnPreconditioner->append( topPelletCladContactPreconditioner );
    }

    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> bottomPelletConstraints;
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> topPelletConstraints;
    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> cladConstraints;
    if ( !useLevitatingFuel ) {
        double fuelOuterRadius = input_db->getDouble( "FuelOuterRadius" );
        double dishRadius      = input_db->getDoubleWithDefault( "DishRadius", 0.0 );
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
        makeConstraintsOnFuel( topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 ),
                               fuelOuterRadius,
                               topPelletConstraints,
                               true,
                               dishRadius );
        makeConstraintsOnFuel( topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 1 ),
                               fuelOuterRadius,
                               topPelletConstraints,
                               false );
    } // end if
    if ( !cladExpansionConstrained ) {
        double cladInnerRadius = input_db->getDouble( "CladInnerRadius" );
        double cladOuterRadius = input_db->getDouble( "CladOuterRadius" );
        double cladHeight      = input_db->getDouble( "CladHeight" );
        makeConstraintsOnClad( cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex ),
                               cladInnerRadius,
                               cladOuterRadius,
                               cladHeight,
                               cladConstraints );
    } // end if

    // Items for computing the RHS correction due to thermal expansion
    AMP::shared_ptr<AMP::Database> fuelTemperatureRhs_db =
        input_db->getDatabase( "FuelTemperatureRHSVectorCorrection" );
    AMP::shared_ptr<AMP::Database> cladTemperatureRhs_db =
        input_db->getDatabase( "CladTemperatureRHSVectorCorrection" );
    AMP::LinearAlgebra::Variable::shared_ptr tempVar(
        new AMP::LinearAlgebra::Variable( "temperature" ) );
    AMP::LinearAlgebra::Variable::shared_ptr dispVar = columnOperator->getOutputVariable();
    AMP::Discretization::DOFManager::shared_ptr tempDofManager =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, 1, split );
    AMP::LinearAlgebra::Vector::shared_ptr tempVec =
        AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    AMP::LinearAlgebra::Vector::shared_ptr refTempVec = tempVec->cloneVector();

    AMP::LinearAlgebra::Vector::shared_ptr sigma_xx = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xx" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_yy = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_yy" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_zz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_zz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_yz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_yz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_xz = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xz" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_xy = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_xy" ) ),
        split );
    AMP::LinearAlgebra::Vector::shared_ptr sigma_eff = AMP::LinearAlgebra::createVector(
        tempDofManager,
        AMP::LinearAlgebra::Variable::shared_ptr( new AMP::LinearAlgebra::Variable( "sigma_eff" ) ),
        split );

    globalComm.barrier();
    double tempCompBeginTime = MPI_Wtime();

    double dummyFuelThermalConductivity = 0.0; // not used k_f = a+b/(c+Td)
    double linearHeatGenerationRate     = input_db->getDouble( "LinearHeatGenerationRate" );
    double fuelOuterRadius              = input_db->getDouble( "FuelOuterRadius" );
    double cladInnerRadius              = input_db->getDouble( "CladInnerRadius" );
    double cladOuterRadius              = input_db->getDouble( "CladOuterRadius" );
    double cladThermalConductivity      = input_db->getDouble( "CladThermalConductivity" );

    double gapThermalConductivity = input_db->getDouble( "GapThermalConductivity" );
    double moderatorTemperature   = input_db->getDouble( "ModeratorTemperature" );
    double moderatorHeatTransferCoefficient =
        input_db->getDouble( "ModeratorHeatTransferCoefficient" );

    double referenceTemperature = input_db->getDouble( "ReferenceTemperature" );

    double cladOuterRadiusTemperature = moderatorTemperature +
                                        linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                            ( cladOuterRadius * moderatorHeatTransferCoefficient );
    double cladInnerRadiusTemperature = cladOuterRadiusTemperature +
                                        linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                            cladThermalConductivity *
                                            std::log( cladOuterRadius / cladInnerRadius );
    double gapAverageRadius = 0.5 * ( cladInnerRadius + fuelOuterRadius );
    double gapHeatTranferCoefficient =
        gapThermalConductivity / ( cladInnerRadius - fuelOuterRadius );
    double fuelOuterRadiusTemperature = cladInnerRadiusTemperature +
                                        linearHeatGenerationRate / ( 2.0 * M_PI ) /
                                            ( gapAverageRadius * gapHeatTranferCoefficient );
    double fuelCenterLineTemperature = fuelOuterRadiusTemperature;
    newton_solver_t<double> solver;
    solver.set( &my_f, &my_ijmf );
    my_p.r       = 0.0;
    my_p.q_prime = linearHeatGenerationRate;
    my_p.T_fo    = fuelOuterRadiusTemperature;
    my_p.R_fo    = fuelOuterRadius;
    solver.solve( fuelCenterLineTemperature, &my_p );
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
            ->getDouble( "THERMAL_EXPANSION_COEFFICIENT" );
    double cladThermalExpansionCoefficient =
        ( cladTemperatureRhs_db->getDatabase( "RhsMaterialModel" ) )
            ->getDouble( "THERMAL_EXPANSION_COEFFICIENT" );
    bottomPelletTopPelletContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );
    bottomPelletCladContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );
    topPelletCladContactOperator->uglyHack(
        tempVec, tempDofManager, fuelThermalExpansionCoefficient, referenceTemperature );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr columnSolVec =
        AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec =
        AMP::LinearAlgebra::createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();
    AMP::LinearAlgebra::Vector::shared_ptr bottomPelletCor;
    AMP::LinearAlgebra::Vector::shared_ptr topPelletCor;
    AMP::LinearAlgebra::Vector::shared_ptr cladCor;

    AMP::LinearAlgebra::Vector::shared_ptr activeSetBeforeUpdateVec = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr activeSetAfterUpdateVec  = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr contactPressureVec       = sigma_eff->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr surfaceTractionVec       = columnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr normalVectorVec          = columnSolVec->cloneVector();

    if ( shrinkFactor != 0.0 ) {
        AMP_ASSERT( ( shrinkFactor > 0.0 ) && ( shrinkFactor < 1.0 ) );
        shrinkMesh( cladMeshAdapter, shrinkFactor );
    }

    if ( fuelNeedALittleHelp != -1.0 ) {
        {
            AMP::Mesh::MeshIterator it =
                topPelletMeshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 2 );
            AMP::Mesh::MeshIterator it_begin = it.begin(), it_end = it.end();
            std::vector<double> coord;
            std::vector<size_t> dofs;
            double epsilon = 1.0e-12;
            double radius;
            for ( it = it_begin; it != it_end; ++it ) {
                coord  = it->coord();
                radius = std::sqrt( std::pow( coord[0], 2 ) + std::pow( coord[1], 2 ) );
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
            AMP::Mesh::MeshIterator it       = cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
            AMP::Mesh::MeshIterator it_begin = it.begin(), it_end = it.end();
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
            AMP::Mesh::MeshIterator it       = cladMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
            AMP::Mesh::MeshIterator it_begin = it.begin(), it_end = it.end();
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

    AMP::LinearAlgebra::Vector::shared_ptr contactShiftVec =
        createVector( dispDofManager, columnVar, split );
    contactShiftVec->zero();

#ifdef USE_EXT_SILO
    {
        siloWriter->registerVector( columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
        siloWriter->registerVector( tempVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Temperature" );
        siloWriter->registerVector( sigma_eff, meshAdapter, AMP::Mesh::GeomType::Vertex, "vonMises" );
        siloWriter->registerVector( sigma_xx, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xx" );
        siloWriter->registerVector( sigma_yy, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_yy" );
        siloWriter->registerVector( sigma_zz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_zz" );
        siloWriter->registerVector( sigma_yz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_yz" );
        siloWriter->registerVector( sigma_xz, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xz" );
        siloWriter->registerVector( sigma_xy, meshAdapter, AMP::Mesh::GeomType::Vertex, "sigma_xy" );
        siloWriter->registerVector(
            activeSetBeforeUpdateVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ActiveSetBeforeUpdate" );
        siloWriter->registerVector(
            activeSetAfterUpdateVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ActiveSetAfterUpdate" );
        siloWriter->registerVector(
            surfaceTractionVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Traction" );
        siloWriter->registerVector( normalVectorVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Normal" );
        siloWriter->registerVector(
            contactPressureVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "ContactPressure" );
        siloWriter->registerVector( contactShiftVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Shift" );
        siloWriter->writeFile( prefixFileName.c_str(), 0 );
    }
#endif

    columnSolVec->zero();
    columnOperator->append( bottomPelletTopPelletContactOperator );
    columnOperator->append( bottomPelletCladContactOperator );
    columnOperator->append( topPelletCladContactOperator );


    // Build a matrix shell operator to use the column operator with the petsc krylov solvers
    AMP::shared_ptr<AMP::Database> matrixShellDatabase =
        input_db->getDatabase( "MatrixShellOperator" );
    AMP::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(
        new AMP::Operator::OperatorParameters( matrixShellDatabase ) );
    AMP::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(
        new AMP::Operator::PetscMatrixShellOperator( matrixShellParams ) );

    int numBottomPelletLocalNodes = 0;
    int numTopPelletLocalNodes    = 0;
    int numCladLocalNodes         = 0;
    if ( bottomPelletMeshAdapter.get() != NULL ) {
        numBottomPelletLocalNodes = bottomPelletMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( topPelletMeshAdapter.get() != NULL ) {
        numTopPelletLocalNodes = topPelletMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    if ( cladMeshAdapter.get() != NULL ) {
        numCladLocalNodes = cladMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
    }
    int matLocalSize =
        dofsPerNode * ( numBottomPelletLocalNodes + numTopPelletLocalNodes + numCladLocalNodes );
    AMP_ASSERT( matLocalSize == static_cast<int>( dispDofManager->numLocalDOF() ) );
    matrixShellOperator->setComm( globalComm );
    matrixShellOperator->setMatLocalRowSize( matLocalSize );
    matrixShellOperator->setMatLocalColumnSize( matLocalSize );
    matrixShellOperator->setOperator( columnOperator );

    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = matrixShellOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = columnPreconditioner;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    AMP::LinearAlgebra::Vector::shared_ptr fullThermalLoadingTempMinusRefTempVec =
        tempVec->cloneVector();
    fullThermalLoadingTempMinusRefTempVec->subtract( tempVec, refTempVec );

    size_t maxActiveSetIterations = input_db->getIntegerWithDefault( "maxActiveSetIterations", 5 );
    size_t maxThermalLoadingIterations =
        input_db->getIntegerWithDefault( "maxThermalLoadingIterations", 5 );
    std::vector<double> scalingFactors;
    std::vector<int> maxIterations;
    bool customLoading = input_db->getBoolWithDefault( "customLoading", false );
    if ( customLoading ) {
        scalingFactors = input_db->getDoubleArray( "scalingFactors" );
        maxIterations  = input_db->getIntegerArray( "maxIterations" );
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
                AMP::LinearAlgebra::Vector::shared_ptr bottomPelletRhsVec =
                    columnRhsVec->select( bottomPelletVectorSelector, dispVar->getName() );
                computeTemperatureRhsVector( bottomPelletMeshAdapter,
                                             fuelTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             bottomPelletRhsVec );

                AMP::LinearAlgebra::VS_Mesh topPelletVectorSelector( topPelletMeshAdapter );
                AMP::LinearAlgebra::Vector::shared_ptr topPelletRhsVec =
                    columnRhsVec->select( topPelletVectorSelector, dispVar->getName() );
                computeTemperatureRhsVector( topPelletMeshAdapter,
                                             fuelTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             topPelletRhsVec );

                AMP::LinearAlgebra::VS_Mesh cladVectorSelector( cladMeshAdapter );
                AMP::LinearAlgebra::Vector::shared_ptr cladRhsVec =
                    columnRhsVec->select( cladVectorSelector, dispVar->getName() );
                computeTemperatureRhsVector( cladMeshAdapter,
                                             cladTemperatureRhs_db,
                                             tempVar,
                                             dispVar,
                                             tempVec,
                                             refTempVec,
                                             cladRhsVec );
            }

            // apply dirichlet rhs correction on f
            if ( bottomPelletBVPOperator.get() != NULL ) {
                bottomPelletBVPOperator->modifyRHSvector( columnRhsVec );
            }
            if ( topPelletBVPOperator.get() != NULL ) {
                topPelletBVPOperator->modifyRHSvector( columnRhsVec );
            }
            if ( cladBVPOperator.get() != NULL ) {
                cladBVPOperator->modifyRHSvector( columnRhsVec );
            }
            {
                AMP::LinearAlgebra::Matrix::shared_ptr bottomPelletMat =
                    bottomPelletBVPOperator->getMatrix();
                AMP::LinearAlgebra::Vector::shared_ptr bottomPelletRhs =
                    bottomPelletBVPOperator->subsetOutputVector( columnRhsVec );
                if ( bottomPelletCor.get() == NULL ) {
                    bottomPelletCor = bottomPelletRhs->cloneVector();
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
                                                   AMP::LinearAlgebra::Matrix::shared_ptr() );
                } // end if
                AMP_ASSERT( bottomPelletCor.get() != NULL );

                AMP::LinearAlgebra::Matrix::shared_ptr topPelletMat =
                    topPelletBVPOperator->getMatrix();
                AMP::LinearAlgebra::Vector::shared_ptr topPelletRhs =
                    topPelletBVPOperator->subsetOutputVector( columnRhsVec );
                if ( topPelletCor.get() == NULL ) {
                    topPelletCor = topPelletRhs->cloneVector();
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
                                                   AMP::LinearAlgebra::Matrix::shared_ptr() );
                } // end if
                AMP_ASSERT( topPelletCor.get() != NULL );

                AMP::LinearAlgebra::Matrix::shared_ptr cladMat = cladBVPOperator->getMatrix();
                AMP::LinearAlgebra::Vector::shared_ptr cladRhs =
                    cladBVPOperator->subsetOutputVector( columnRhsVec );
                if ( cladCor.get() == NULL ) {
                    cladCor = cladRhs->cloneVector();
                    applyCustomDirichletCondition(
                        cladRhs, cladCor, cladMeshAdapter, cladConstraints, cladMat );
                } else {
                    applyCustomDirichletCondition( cladRhs,
                                                   cladCor,
                                                   cladMeshAdapter,
                                                   cladConstraints,
                                                   AMP::LinearAlgebra::Matrix::shared_ptr() );
                } // end if
                AMP_ASSERT( topPelletCor.get() != NULL );
            }

            // get d
            contactShiftVec->zero();
            bottomPelletTopPelletContactOperator->addShiftToSlave( contactShiftVec );
            bottomPelletCladContactOperator->addShiftToSlave( contactShiftVec );
            topPelletCladContactOperator->addShiftToSlave( contactShiftVec );

            // compute - Kd
            AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec =
                createVector( dispDofManager, columnVar, split );
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

            linearSolver->solve( columnRhsVec, columnSolVec );

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

            std::vector<AMP::Mesh::MeshElementID> const &bottomPelletTopPelletActiveSet =
                bottomPelletTopPelletContactOperator->getActiveSet();
            size_t const bottomPelletTopPelletSizeOfActiveSetBeforeUpdate =
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
            std::vector<AMP::Mesh::MeshElementID> const &bottomPelletCladActiveSet =
                bottomPelletCladContactOperator->getActiveSet();
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
            std::vector<AMP::Mesh::MeshElementID> const &topPelletCladActiveSet =
                topPelletCladContactOperator->getActiveSet();
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


#ifdef USE_EXT_SILO
            {

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
            }
#endif

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


#ifdef USE_EXT_SILO
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
#endif

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
