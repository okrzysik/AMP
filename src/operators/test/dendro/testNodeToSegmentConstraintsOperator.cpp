#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshParameters.h"
#include "AMP/ampmesh/euclidean_geometry_tools.h"
#include "AMP/ampmesh/latex_visualization_tools.h"
#include "AMP/ampmesh/libmesh/libmeshMesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/contact/NodeToSegmentConstraintsOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"
#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/contact/MPCSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/ReadTestMesh.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "externVars.h"

#include <fstream>


static void drawVerticesOnBoundaryID( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                      int boundaryID,
                                      std::ostream &os,
                                      double const *point_of_view,
                                      const std::string &option = "" )
{
    auto boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID );
    auto boundaryIterator_begin = boundaryIterator.begin();
    auto boundaryIterator_end   = boundaryIterator.end();
    std::vector<double> vertexCoordinates;
    os << std::setprecision( 6 ) << std::fixed;
    for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
          ++boundaryIterator ) {
        vertexCoordinates = boundaryIterator->coord();
        AMP_ASSERT( vertexCoordinates.size() == 3 );
        draw_point( &( vertexCoordinates[0] ), option, os );
    } // end for
}

static void drawGeomType::FacesOnBoundaryID( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                             int boundaryID,
                                             std::ostream &os,
                                             double const *point_of_view,
                                             const std::string &option = "" )
{
    auto boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, boundaryID );
    auto boundaryIterator_begin = boundaryIterator.begin();
    auto boundaryIterator_end   = boundaryIterator.end();
    std::vector<AMP::Mesh::MeshElement> faceVertices;
    std::vector<double> faceGeomType::VertexCoordinates;
    double faceData[12]          = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double const *faceDataPtr[4] = { faceData, faceData + 3, faceData + 6, faceData + 9 };
    os << std::setprecision( 6 ) << std::fixed;
    for ( boundaryIterator = boundaryIterator_begin; boundaryIterator != boundaryIterator_end;
          ++boundaryIterator ) {
        faceVertices = boundaryIterator->getElements( AMP::Mesh::GeomType::Vertex );
        AMP_ASSERT( faceVertices.size() == 4 );
        for ( size_t i = 0; i < 4; ++i ) {
            faceGeomType::VertexCoordinates = faceVertices[i].coord();
            AMP_ASSERT( faceGeomType::VertexCoordinates.size() == 3 );
            std::copy( faceGeomType::VertexCoordinates.begin(),
                       faceGeomType::VertexCoordinates.end(),
                       faceData + 3 * i );
        } // end for i
        triangle_t t( faceDataPtr[0], faceDataPtr[1], faceDataPtr[2] );

        if ( compute_scalar_product( point_of_view, t.get_normal() ) > 0.0 ) {
            os << "\\draw[" << option << "]\n";
            write_face( faceDataPtr, os );
        } // end if
    }     // end for
}

static void myPCG( AMP::LinearAlgebra::Vector::shared_ptr rhs,
                   AMP::LinearAlgebra::Vector::shared_ptr sol,
                   AMP::Operator::Operator::shared_ptr op,
                   std::shared_ptr<AMP::Solver::SolverStrategy> pre,
                   size_t maxIters,
                   double relTol,
                   double absTol,
                   bool verbose     = false,
                   std::ostream &os = std::cout )
{
    auto res    = sol->cloneVector();
    auto dir    = sol->cloneVector();
    auto ext    = sol->cloneVector();
    auto oldSol = sol->cloneVector();
    auto oldRes = sol->cloneVector();
    auto oldDir = sol->cloneVector();
    auto matVec = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    op->apply( nullVec, sol, matVec, 1.0, 0.0 );
    oldRes->subtract( rhs, matVec );
    pre->solve( oldRes, ext );
    oldDir->copyVector( ext );
    oldSol->copyVector( sol );
    double initialResNorm = oldRes->L2Norm();
    double tol            = absTol + relTol * initialResNorm;
    if ( verbose ) {
        os << std::setprecision( 15 ) << "  iter=0  itialResNorm=" << initialResNorm << "\n";
    }
    for ( size_t iter = 0; iter < maxIters; ++iter ) {
        if ( verbose ) {
            os << "  iter=" << iter + 1 << "  ";
        }
        op->apply( nullVec, oldDir, matVec, 1.0, 0.0 );
        double extDOToldRes    = ext->dot( oldRes );
        double oldDirDOTmatVec = oldDir->dot( matVec );
        double alpha           = extDOToldRes / oldDirDOTmatVec;
        if ( verbose ) {
            os << "alpha=" << alpha << "  ";
        }
        if ( verbose ) {
            os << "oldDirDOTmatVec=" << oldDirDOTmatVec << "  ";
        }
        sol->axpy( alpha, oldDir, oldSol );
        res->axpy( -alpha, matVec, oldRes );
        double resNorm = res->L2Norm();
        if ( verbose ) {
            os << "resNorm=" << resNorm << "  ";
        }
        if ( resNorm < tol ) {
            os << "\n";
            break;
        }
        pre->solve( res, ext );
        double extDOTres = ext->dot( res );
        double beta      = extDOTres / extDOToldRes;
        if ( verbose ) {
            os << "beta=" << beta << "\n";
        }
        dir->axpy( beta, oldDir, ext );
        oldSol->copyVector( sol );
        oldRes->copyVector( res );
        oldDir->copyVector( dir );
    } // end for
}


static void myTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
#endif

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    std::fstream fout;
    std::string filename = "debug_driver_" + amp::utilities::inttostring( rank );
    fout.open( filename.c_str(), std::fstream::out );

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
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

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
    auto dofManager     = AMP::Discretization::simpleDOFManager::create(
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

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    auto contact_db = input_db->getDatabase( "ContactOperator" );
    auto contactOperatorParams =
        std::make_shared<AMP::Operator::NodeToSegmentConstraintsOperatorParameters>( contact_db );
    contactOperatorParams->d_DOFsPerNode = dofsPerNode;
    contactOperatorParams->d_DOFManager  = dofManager;
    contactOperatorParams->d_GlobalComm  = globalComm;
    contactOperatorParams->d_Mesh        = meshAdapter;

    auto contactOperator =
        std::make_shared<AMP::Operator::NodeToSegmentConstraintsOperator>( contactOperatorParams );

    // TODO: RESET IN CONSTRUCTOR?
    contactOperator->reset( contactOperatorParams );

    // Build the master and slave operators
    auto masterMeshID      = contactOperator->getMasterMeshID();
    auto masterMeshAdapter = meshAdapter->Subset( masterMeshID );
    if ( masterMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
        auto masterOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
        columnOperator->append( masterOperator );

        auto masterSolver_db = columnPreconditioner_db->getDatabase( "MasterSolver" );
        auto masterSolverParams =
            std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( masterSolver_db );
        masterSolverParams->d_pOperator = masterOperator;
        masterSolverParams->d_comm      = masterMeshAdapter->getComm();
        //    masterSolverParams->d_comm = globalComm;
        auto masterSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( masterSolverParams );
        columnPreconditioner->append( masterSolver );
        std::fstream masterFout;
        masterFout.open( "master_pellet", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 1, masterFout, point_of_view );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 4, masterFout, point_of_view );
        masterFout.close();
    } // end if

    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

    auto slaveMeshID      = contactOperator->getSlaveMeshID();
    auto slaveMeshAdapter = meshAdapter->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        std::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;

        auto slaveSolver_db = columnPreconditioner_db->getDatabase( "SlaveSolver" );
        auto slaveSolverParams =
            std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( slaveSolver_db );

        bool useSlaveBVPOperator = input_db->getScalar<bool>( "useSlaveBVPOperator" );
        if ( useSlaveBVPOperator ) {
            slaveBVPOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel ) );
            columnOperator->append( slaveBVPOperator );
            slaveSolverParams->d_pOperator = slaveBVPOperator;
        } else {
            auto slaveMechanicsLinearFEOperator =
                std::dynamic_pointer_cast<AMP::Operator::MechanicsLinearFEOperator>(
                    AMP::Operator::OperatorBuilder::createOperator( slaveMeshAdapter,
                                                                    "MechanicsLinearFEOperator",
                                                                    input_db,
                                                                    slaveElementPhysicsModel ) );
            columnOperator->append( slaveMechanicsLinearFEOperator );

            slaveSolverParams->d_pOperator = slaveMechanicsLinearFEOperator;

            slaveLoadOperator = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    slaveMeshAdapter, "SlaveLoadOperator", input_db, slaveElementPhysicsModel ) );
            auto slaveVar = slaveMechanicsLinearFEOperator->getOutputVariable();
            slaveLoadOperator->setVariable( slaveVar );
        } // end if

        //    slaveSolverParams->d_comm = globalComm;
        slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
        auto slaveSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( slaveSolverParams );
        columnPreconditioner->append( slaveSolver );

        std::fstream slaveFout;
        slaveFout.open( "slave_pellet", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed" );
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed" );
        drawVerticesOnBoundaryID( slaveMeshAdapter, 2, slaveFout, point_of_view, "red" );
        slaveFout.close();
    } // end if

    auto contactPreconditioner_db = columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    auto contactPreconditionerParams =
        std::make_shared<AMP::Solver::MPCSolverParameters>( contactPreconditioner_db );
    contactPreconditionerParams->d_pOperator = contactOperator;
    auto contactPreconditioner =
        std::make_shared<AMP::Solver::MPCSolver>( contactPreconditionerParams );
    columnPreconditioner->append( contactPreconditioner );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnVar    = columnOperator->getOutputVariable();
    auto columnSolVec = createVector( dofManager, columnVar, split );
    auto columnRhsVec = createVector( dofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();

    // compute f
    if ( slaveLoadOperator.get() != NULL ) {
        slaveLoadOperator->apply( nullVec, nullVec, columnRhsVec, 1.0, 0.0 );
    } // end if

    // apply dirichlet rhs correction
    if ( slaveBVPOperator.get() != NULL ) {
        slaveBVPOperator->modifyRHSvector( columnRhsVec );
    } // end if

    // get d
    contactOperator->addShiftToSlave( columnSolVec );

    // compute - Kd
    auto rhsCorrectionVec = createVector( dofManager, columnVar, split );
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

    bool usePetscKrylovSolver = input_db->getScalar<bool>( "usePetscKrylovSolver" );
    if ( usePetscKrylovSolver ) {
        // Build a matrix shell operator to use the column operator with the petsc krylov solvers
        auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
        auto matrixShellParams =
            std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase );
        auto matrixShellOperator =
            std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );

        int numMasterLocalNodes = 0;
        int numSlaveLocalNodes  = 0;
        if ( masterMeshAdapter.get() != NULL ) {
            numMasterLocalNodes =
                masterMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
        }
        if ( slaveMeshAdapter.get() != NULL ) {
            numSlaveLocalNodes = slaveMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
        }
        int matLocalSize = dofsPerNode * ( numMasterLocalNodes + numSlaveLocalNodes );
        AMP_ASSERT( matLocalSize == static_cast<int>( dofManager->numLocalDOF() ) );
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

        linearSolver->solve( columnRhsVec, columnSolVec );
    } else {
        size_t myPCGmaxIters = input_db->getScalar<int>( "myPCGmaxIters" );
        double myPCGrelTol   = input_db->getScalar<double>( "myPCGrelTol" );
        double myPCGabsTol   = input_db->getScalar<double>( "myPCGabsTol" );
        myPCG( columnRhsVec,
               columnSolVec,
               columnOperator,
               columnPreconditioner,
               myPCGmaxIters,
               myPCGrelTol,
               myPCGabsTol,
               true );
    }
    // u^s = C u^m + d
    contactOperator->copyMasterToSlave( columnSolVec );
    contactOperator->addShiftToSlave( columnSolVec );

    meshAdapter->displaceMesh( columnSolVec );

    if ( masterMeshAdapter.get() != NULL ) {
        std::fstream masterFout;
        masterFout.open( "master_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 1, masterFout, point_of_view, "" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 4, masterFout, point_of_view, "" );
        masterFout.close();
    } // end if
    if ( slaveMeshAdapter.get() != NULL ) {
        std::fstream slaveFout;
        slaveFout.open( "slave_pellet_displaced_mesh", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed" );
        drawGeomType::FacesOnBoundaryID( slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed" );
        // drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        slaveFout.close();
    } // end if

#ifdef USE_EXT_SILO
    siloWriter->registerVector(
        columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
    char outFileName[256];
    sprintf( outFileName, "MPC_%d", 0 );
    siloWriter->writeFile( outFileName, 0 );
#endif
    fout.close();

    ut->passes( exeName );
}

static void myTest2( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    //  int npes = globalComm.getSize();
    int rank = globalComm.getRank();
    if ( !rank ) {
        std::cout << "### FUSED MESHES CASE ###" << std::endl;
    }

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

    bool skipFusedMesh = input_db->getScalar<bool>( "skipFusedMesh" );
    if ( skipFusedMesh ) {
        return;
    }

    // Load the meshes
    globalComm.barrier();
    double meshBeginTime = MPI_Wtime();

    AMP_INSIST( input_db->keyExists( "FusedMesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db    = input_db->getDatabase( "FusedMesh" );
    auto meshParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    meshParams->setComm( globalComm );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( meshParams );

    globalComm.barrier();
    double meshEndTime = MPI_Wtime();
    if ( !rank ) {
        std::cout << "Finished reading the mesh in " << ( meshEndTime - meshBeginTime )
                  << " seconds." << std::endl;
    }


    int dofsPerNode     = 3;
    int nodalGhostWidth = 1;
    bool split          = true;
    auto dofManager     = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, dofsPerNode, split );


    // build a column operator and a column preconditioner
    auto columnOperator          = std::make_shared<AMP::Operator::ColumnOperator>();
    auto linearSolver_db         = input_db->getDatabase( "LinearSolver" );
    auto columnPreconditioner_db = linearSolver_db->getDatabase( "Preconditioner" );
    auto columnPreconditionerParams =
        std::make_shared<AMP::Solver::ColumnSolverParameters>( columnPreconditioner_db );
    columnPreconditionerParams->d_pOperator = columnOperator;
    auto columnPreconditioner =
        std::make_shared<AMP::Solver::ColumnSolver>( columnPreconditionerParams );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> masterOperator;
    bool useSlaveBVPOperator = input_db->getScalar<bool>( "useSlaveBVPOperator" );
    if ( useSlaveBVPOperator ) {
        masterOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "FusedMeshBVPOperator", input_db, masterElementPhysicsModel ) );
    } else {
        masterOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
    }
    columnOperator->append( masterOperator );

    auto masterSolver_db = columnPreconditioner_db->getDatabase( "MasterSolver" );
    auto masterSolverParams =
        std::make_shared<AMP::Solver::PetscKrylovSolverParameters>( masterSolver_db );
    masterSolverParams->d_pOperator = masterOperator;
    masterSolverParams->d_comm      = globalComm;
    auto masterSolver = std::make_shared<AMP::Solver::PetscKrylovSolver>( masterSolverParams );
    columnPreconditioner->append( masterSolver );

    auto slaveLoadOperator = std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "SlaveLoadOperator", input_db, masterElementPhysicsModel ) );
    auto slaveVar = masterOperator->getOutputVariable();
    slaveLoadOperator->setVariable( slaveVar );


    auto columnVar = columnOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    auto columnSolVec = createVector( dofManager, columnVar, split );
    auto columnRhsVec = createVector( dofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();

    if ( useSlaveBVPOperator ) {
        masterOperator->modifyRHSvector( columnRhsVec );
    } else {
        slaveLoadOperator->apply( nullVec, nullVec, columnRhsVec, 1.0, 0.0 );
    }

    bool usePetscKrylovSolver = input_db->getScalar<bool>( "usePetscKrylovSolver" );
    if ( usePetscKrylovSolver ) {
        // Build a matrix shell operator to use the column operator with the petsc krylov solvers
        auto matrixShellDatabase = input_db->getDatabase( "MatrixShellOperator" );
        auto matrixShellParams =
            std::make_shared<AMP::Operator::OperatorParameters>( matrixShellDatabase );
        auto matrixShellOperator =
            std::make_shared<AMP::Operator::PetscMatrixShellOperator>( matrixShellParams );

        int matLocalSize =
            dofsPerNode * meshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
        AMP_ASSERT( matLocalSize == static_cast<int>( dofManager->numLocalDOF() ) );
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
        linearSolver->setZeroInitialGuess( true );

        linearSolver->solve( columnRhsVec, columnSolVec );
    } else {
        size_t myPCGmaxIters = input_db->getScalar<int>( "myPCGmaxIters" );
        double myPCGrelTol   = input_db->getScalar<double>( "myPCGrelTol" );
        double myPCGabsTol   = input_db->getScalar<double>( "myPCGabsTol" );
        myPCG( columnRhsVec,
               columnSolVec,
               columnOperator,
               columnPreconditioner,
               myPCGmaxIters,
               myPCGrelTol,
               myPCGabsTol,
               true );
    }

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    //  exeNames.push_back("testNodeToSegmentConstraintsOperator-cube");
    //  exeNames.push_back("testNodeToSegmentConstraintsOperator-cylinder");
    exeNames.push_back( "testNodeToSegmentConstraintsOperator-pellet" );

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
        myTest2( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
