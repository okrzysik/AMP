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
#include "AMP/operators/CustomConstraintsEliminationOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/TrilinosMatrixShellOperator.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/contact/NodeToGeomType::FaceContactOperator.h"
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


static void selectNodes( AMP::Mesh::Mesh::shared_ptr mesh,
                         std::vector<AMP::Mesh::MeshElementID> &nodesGlobalIDs )
{
    auto meshIterator       = mesh->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 3 );
    auto meshIterator_begin = meshIterator.begin();
    auto meshIterator_end   = meshIterator.end();
    nodesGlobalIDs.clear();
    for ( meshIterator = meshIterator_begin; meshIterator != meshIterator_end; ++meshIterator ) {
        std::vector<double> coord = meshIterator->coord();
        if ( std::abs( coord[2] - 0.0 ) < 1.0e-14 ) {
            nodesGlobalIDs.push_back( meshIterator->globalID() );
            //      std::cout<<nodesGlobalIDs.size()<<"  ("<<coord[0]<<", "<<coord[1]<<",
            //      "<<coord[2]<<")"<<std::endl;
        } // end if
    }     // end for
}

static void printNodesValues( AMP::Mesh::Mesh::shared_ptr mesh,
                              std::vector<AMP::Mesh::MeshElementID> const &nodesGlobalIDs,
                              AMP::LinearAlgebra::Vector::shared_ptr vectorField,
                              std::ostream &os = std::cout )
{
    auto dofManager = vectorField->getDOFManager();
    for ( size_t i = 0; i < nodesGlobalIDs.size(); ++i ) {
        auto coord = mesh->getElement( nodesGlobalIDs[i] ).coord();
        std::vector<size_t> dofIndices;
        dofManager->getDOFs( nodesGlobalIDs[i], dofIndices );
        AMP_ASSERT( dofIndices.size() == 1 );
        double value = vectorField->getLocalValueByGlobalID( dofIndices[0] );
        os << std::setprecision( 15 ) << coord[0] << "  " << value << "\n";
    } // end for i
}

static void
getConcentratedLoadAtNodes( double loadParameter,
                            double loadCutoff,
                            AMP::Mesh::Mesh::shared_ptr meshAdapter,
                            AMP::LinearAlgebra::Vector::shared_ptr loadVector,
                            std::shared_ptr<AMP::Discretization::DOFManager> dofManager )
{
    static std::vector<double> loadValues;
    static std::vector<size_t> dofIndices;
    AMP_ASSERT( loadValues.size() == dofIndices.size() );

    if ( loadValues.empty() ) {
        double totalLoad = 0.0;
        auto boundaryIterator =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, 4, 0 );
        auto boundaryIterator_begin = boundaryIterator.begin();
        auto boundaryIterator_end   = boundaryIterator.end();
        std::vector<double> vertexCoordinates;
        std::vector<size_t> vertexDofIndices;
        for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
              ++boundaryIterator ) {
            vertexCoordinates = boundaryIterator->coord();
            AMP_ASSERT( vertexCoordinates.size() == 3 );
            dofManager->getDOFs( boundaryIterator->globalID(), vertexDofIndices );
            AMP_ASSERT( vertexDofIndices.size() == 3 );

            if ( vertexCoordinates[1] > loadCutoff ) {
                loadValues.push_back( loadParameter );
                dofIndices.push_back( vertexDofIndices[1] );
                if ( vertexCoordinates[1] - 1.98426 < 0.0 ) {
                    loadValues.back() /= 2.0;
                } // end if
                if ( ( std::abs( vertexCoordinates[2] + 1.0 ) < 1.0e-14 ) ||
                     ( std::abs( vertexCoordinates[2] - 1.0 ) < 1.0e-14 ) ) {
                    loadValues.back() /= 2.0;
                } // end if
                totalLoad += loadValues.back();
                //        std::cout<<loadValues.size()<<"  "<<loadValues.back()<<"
                //        "<<dofIndices.back()<<"
                //        ("<<vertexCoordinates[0]<<", "<<vertexCoordinates[1]<<"
                //        ,"<<vertexCoordinates[2]<<")\n";
            } // end if
        }     // end for
        std::cout << "TOTAL load=" << totalLoad << "\n";
        AMP_ASSERT( loadValues.size() > 0 );
    } // end if

    loadVector->zero();
    loadVector->setLocalValuesByGlobalID(
        loadValues.size(), &( dofIndices[0] ), &( loadValues[0] ) );
    loadVector->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
}

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

    // Get the mechanics material model for the contact operator
    auto model_db = input_db->getDatabase( "MasterMechanicsMaterialModel" );
    auto masterMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( model_db );
    auto masterMechanicsMaterialModel = std::make_shared<AMP::Operator::IsotropicElasticModel>(
        masterMechanicsMaterialModelParams );

    // ... needed for computing stresses
    auto slaveMechanicsMaterialModel_db = input_db->getDatabase( "SlaveMechanicsMaterialModel" );
    auto slaveMechanicsMaterialModelParams =
        std::make_shared<AMP::Operator::MechanicsModelParameters>( slaveMechanicsMaterialModel_db );
    auto slaveMechanicsMaterialModel =
        std::make_shared<AMP::Operator::IsotropicElasticModel>( slaveMechanicsMaterialModelParams );

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    auto contact_db = input_db->getDatabase( "ContactOperator" );
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
        contact_db->getWithDefault<bool>( "ContactIsFrictionless", false ) );

    bool useML = input_db->getWithDefault<bool>( "useML", false );

    // Build the master and slave operators
    std::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> dummyMasterBVPOperator;
    auto masterMeshID      = contactOperator->getMasterMeshID();
    auto masterMeshAdapter = meshAdapter->Subset( masterMeshID );
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

    } // end if

    std::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;
    auto slaveMeshID      = contactOperator->getSlaveMeshID();
    auto slaveMeshAdapter = meshAdapter->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;
        slaveBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel ) );
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

    } // end if

    std::map<AMP::Mesh::MeshElementID, std::map<size_t, double>> constraints;
    auto it = masterMeshAdapter->getIterator( AMP::Mesh::GeomType::Vertex );
    {
        AMP::Mesh::MeshIterator it_begin = it.begin();
        AMP::Mesh::MeshIterator it_end   = it.end();
        std::vector<double> coord;
        double const epsilon = 1.0e-10;
        std::map<size_t, double> tmp;
        for ( it = it_begin; it != it_end; ++it ) {
            coord = it->coord();
            tmp.clear();
            if ( std::abs(
                     std::sqrt( std::pow( coord[0] - 0.0, 2 ) + std::pow( coord[1] - 1.0, 2 ) ) -
                     1.0 ) < epsilon ) {
                if ( ( std::abs( coord[2] - 0.0 ) < epsilon ) && ( coord[1] > 1.0 - epsilon ) ) {
                    tmp.insert( std::pair<size_t, double>( 2, 0.0 ) );
                }
                if ( std::abs( coord[1] - 0.0 ) < epsilon ) {
                    tmp.insert( std::pair<size_t, double>( 1, 0.0 ) );
                }
            } else {
                if ( ( ( std::abs( coord[2] + 1.0 ) < epsilon ) ||
                       ( std::abs( coord[2] - 1.0 ) < epsilon ) ) &&
                     ( std::abs( coord[0] - 0.0 ) < epsilon ) ) {
                    tmp.insert( std::pair<size_t, double>( 0, 0.0 ) );
                }
            } // end if
            if ( !tmp.empty() ) {
                constraints.insert( std::pair<AMP::Mesh::MeshElementID, std::map<size_t, double>>(
                    it->globalID(), tmp ) );
            } // end if
        }     // end for
    }

    // Build matrix shell operators to use the column operator with the petsc krylov solvers and
    // trilinos ml
    auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
    auto matrixShellParams =
        std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase );

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

    auto petscMatrixShellOperator =
        std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );
    petscMatrixShellOperator->setComm( globalComm );
    petscMatrixShellOperator->setMatLocalRowSize( matLocalSize );
    petscMatrixShellOperator->setMatLocalColumnSize( matLocalSize );
    petscMatrixShellOperator->setOperator( columnOperator );

    auto contactPreconditioner_db = columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    auto contactPreconditionerParams =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolverParameters>(
            contactPreconditioner_db );
    contactPreconditionerParams->d_pOperator = contactOperator;
    auto contactPreconditioner =
        std::make_shared<AMP::Solver::ConstraintsEliminationSolver>( contactPreconditionerParams );
    columnPreconditioner->append( contactPreconditioner );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnVar    = columnOperator->getOutputVariable();
    auto columnSolVec = createVector( dispDofManager, columnVar, split );
    auto columnRhsVec = createVector( dispDofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();
    AMP::LinearAlgebra::Vector::shared_ptr cor;
    AMP_ASSERT( cor.get() == nullptr );

    auto tempVar        = std::make_shared<MP::LinearAlgebra::Variable>( "temperature" );
    auto dispVar        = columnOperator->getOutputVariable();
    auto tempDofManager = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, 1, split );
    auto tempVec = AMP::LinearAlgebra::createVector( tempDofManager, tempVar, split );
    double const referenceTemperature = 300.0;
    tempVec->setToScalar( referenceTemperature );
    double const thermalExpansionCoefficient = 2.0e-4;

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
    auto activeSetBeforeUpdateVec = sigma_eff->cloneVector();
    auto activeSetAfterUpdateVec  = sigma_eff->cloneVector();
    auto contactPressureVec       = sigma_eff->cloneVector();
    auto surfaceTractionVec       = columnSolVec->cloneVector();
    auto normalVectorVec          = columnSolVec->cloneVector();
    auto contactShiftVec          = columnSolVec->cloneVector();
    contactPressureVec->zero();
    surfaceTractionVec->zero();
    normalVectorVec->zero();
    contactShiftVec->zero();
    sigma_xx->zero();
    sigma_yy->zero();
    sigma_zz->zero();
    sigma_yz->zero();
    sigma_xz->zero();
    sigma_xy->zero();

    bool skipDisplaceMesh = true;
    contactOperator->updateActiveSet( nullVec, skipDisplaceMesh );

    siloWriter->registerVector(
        columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "SolutionDisplacement" );
    siloWriter->registerVector(
        sigma_eff, meshAdapter, AMP::Mesh::GeomType::Vertex, "vonMisesStresses" );
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
    siloWriter->writeFile( "TITI_0", 0 );

    auto linearSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( linearSolver_db );
    linearSolverParams->d_pOperator       = petscMatrixShellOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = columnPreconditioner;
    auto linearSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( linearSolverParams );
    //  linearSolver->setZeroInitialGuess(true);
    linearSolver->setInitialGuess( columnSolVec );

    //  std::vector<AMP::Mesh::MeshElementID> masterNodesGlobalIDs;
    //  selectNodes(masterMeshAdapter, masterNodesGlobalIDs);
    //  printNodesValues(masterMeshAdapter, masterNodesGlobalIDs, contactPressureVec);
    std::vector<AMP::Mesh::MeshElementID> slaveNodesGlobalIDs;
    selectNodes( slaveMeshAdapter, slaveNodesGlobalIDs );
    printNodesValues( slaveMeshAdapter, slaveNodesGlobalIDs, contactPressureVec );

    int TOTO_count = 0;
    size_t const maxLoadingIterations =
        input_db->getWithDefault<size_t>( "maxLoadingIterations", 5 );
    for ( size_t loadingIteration = 0; loadingIteration < maxLoadingIterations;
          ++loadingIteration ) {
        double scalingFactor = static_cast<double>( loadingIteration + 1 ) /
                               static_cast<double>( maxLoadingIterations );
        // quadratic scheme
        //  scalingFactor = scalingFactor * scalingFactor;
        if ( !rank ) {
            std::cout << "LOADING STEP " << loadingIteration + 1 << "/" << maxLoadingIterations
                      << "  (factor=" << scalingFactor << ")\n";
        }

        size_t const maxActiveSetIterations =
            input_db->getWithDefault<size_t>( "maxActiveSetIterations", 5 );
        for ( size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations;
              ++activeSetIteration ) {
            if ( !rank ) {
                std::cout << "ACTIVE SET ITERATION #" << activeSetIteration + 1 << "\n";
            }
            ++TOTO_count;

            columnSolVec->zero();
            columnRhsVec->zero();

            // compute rhs
            double loadParameter = input_db->getScalar<double>( "LoadParameter" );
            double loadCutoff    = input_db->getScalar<double>( "LoadCutoff" );
            //  getConcentratedLoadAtNodes(loadParameter, loadCutoff, slaveMeshAdapter,
            //  columnRhsVec, dispDofManager);
            getConcentratedLoadAtNodes(
                loadParameter, loadCutoff, masterMeshAdapter, columnRhsVec, dispDofManager );
            columnRhsVec->scale( scalingFactor );

            // apply dirichlet rhs correction
            if ( masterBVPOperator.get() != NULL ) {
                masterBVPOperator->modifyRHSvector( columnRhsVec );
            } // end if
            if ( slaveBVPOperator.get() != NULL ) {
                slaveBVPOperator->modifyRHSvector( columnRhsVec );
            } // end if

            auto mat = masterBVPOperator->getMatrix();
            auto rhs = masterBVPOperator->subsetOutputVector( columnRhsVec );
            if ( cor.get() == nullptr ) {
                cor = rhs->cloneVector();
                applyCustomDirichletCondition( rhs, cor, meshAdapter, constraints, mat );
            } else {
                applyCustomDirichletCondition( rhs,
                                               cor,
                                               meshAdapter,
                                               constraints,
                                               std::shared_ptr<AMP::LinearAlgebra::Matrix>() );
            } // end if
            AMP_ASSERT( cor.get() != nullptr );

            // get d
            contactOperator->addShiftToSlave( columnSolVec );

            // compute - Kd
            auto rhsCorrectionVec = createVector( dispDofManager, columnVar, split );
            columnOperator->apply( nullVec, columnSolVec, rhsCorrectionVec, -1.0, 0.0 );
            columnOperator->append( contactOperator );

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
                                 thermalExpansionCoefficient,
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
                                 thermalExpansionCoefficient,
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

            // Update active set
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

            //{
            //  auto meshIterator =
            //  masterMeshAdapter->getBoundaryIDIterator(AMP::Mesh::GeomType::Vertex, 3);
            //  auto meshIterator_begin = meshIterator.begin();
            //  auto meshIterator_end = meshIterator.end();
            //  double normalVector[3] = { 0.0, 1.0, 0.0 };
            //  double stressTensor[6];
            //  double traction[3];
            //  double contactPressure;
            //  std::vector<size_t> dofIndex;
            //  for (meshIterator = meshIterator_begin; meshIterator != meshIterator_end;
            //  ++meshIterator) {
            //    tempDofManager->getDOFs(meshIterator->globalID(), dofIndex);
            //    AMP_ASSERT(dofIndex.size() == 1);
            //    stressTensor[0] = sigma_xx->getLocalValueByGlobalID(dofIndex[0]);
            //    stressTensor[1] = sigma_yy->getLocalValueByGlobalID(dofIndex[0]);
            //    stressTensor[2] = sigma_zz->getLocalValueByGlobalID(dofIndex[0]);
            //    stressTensor[3] = sigma_yz->getLocalValueByGlobalID(dofIndex[0]);
            //    stressTensor[4] = sigma_xz->getLocalValueByGlobalID(dofIndex[0]);
            //    stressTensor[5] = sigma_xy->getLocalValueByGlobalID(dofIndex[0]);
            //    compute_traction(stressTensor, normalVector, traction);
            //    contactPressure = - compute_scalar_product(traction, normalVector);
            //    contactPressureVec->setLocalValueByGlobalID(dofIndex[0], contactPressure);
            //  } // end for
            //}

            //  printNodesValues(masterMeshAdapter, masterNodesGlobalIDs, contactPressureVec);
            printNodesValues( slaveMeshAdapter, slaveNodesGlobalIDs, contactPressureVec );

            //  std::vector<double> coord;
            //  std::vector<size_t> dof;
            //  double disp[3];
            //  for (size_t jj = 0; jj < sizeOfActiveSetBeforeUpdate; ++jj) {
            //    auto coord = (slaveMeshAdapter->getElement(activeSet[jj])).coord();
            //    if (std::abs(coord[2] - 0.0) < 1.0e-10) {
            //      tempDofManager->getDOFs(activeSet[jj], dof);
            //      double pressure = contactPressureVec->getLocalValueByGlobalID(dof[0]);
            //      dispDofManager->getDOFs(activeSet[jj], dof);
            //      columnSolVec->getLocalValuesByGlobalID(3, &(dof[0]), disp);
            //      std::transform(&(coord[0]), &(coord[0])+3, &(coord[0]), disp,
            //      std::plus<double>()); std::cout<<coord[0]<<"  "<<pressure<<"\n";
            //    } // end if
            //  } // end for jj

            meshAdapter->displaceMesh( columnSolVec );
            siloWriter->writeFile( "TITI_0", TOTO_count );
            columnSolVec->scale( -1.0 );
            meshAdapter->displaceMesh( columnSolVec );
            columnSolVec->scale( -1.0 );

            if ( !rank )
                std::cout << nChangesInActiveSet << " CHANGES IN ACTIVE SET\n";

            if ( nChangesInActiveSet == 0 )
                break;
            AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
        } // end for

    } // end for

    fout.close();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNodeToGeomType::FaceContactOperator-3" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
