
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/InputDatabase.h"
#include "AMP/utils/InputManager.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include "externVars.h"

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/libmesh/libMesh.h"
#include "AMP/utils/Writer.h"

#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/contact/NodeToGeomType::FaceContactOperator.h"
#include "AMP/operators/mechanics/MechanicsLinearFEOperator.h"
#include "AMP/operators/mechanics/MechanicsMaterialModel.h"
#include "AMP/operators/mechanics/MechanicsModelParameters.h"
#include "AMP/operators/petsc/PetscMatrixShellOperator.h"

#include "AMP/solvers/ColumnSolver.h"
#include "AMP/solvers/ConstraintsEliminationSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscKrylovSolverParameters.h"

#include "AMP/utils/ReadTestMesh.h"

#include "AMP/ampmesh/euclidean_geometry_tools.h"
#include "AMP/ampmesh/latex_visualization_tools.h"
#include <fstream>


static void drawVerticesOnBoundaryID( AMP::Mesh::Mesh::shared_ptr meshAdapter,
                                      int boundaryID,
                                      std::ostream &os,
                                      double const *point_of_view,
                                      const std::string &option = "" )
{
    AMP::Mesh::MeshIterator boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, boundaryID );
    AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(),
                            boundaryIterator_end   = boundaryIterator.end();
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
    AMP::Mesh::MeshIterator boundaryIterator =
        meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Face, boundaryID );
    AMP::Mesh::MeshIterator boundaryIterator_begin = boundaryIterator.begin(),
                            boundaryIterator_end   = boundaryIterator.end();
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

        //    if (compute_scalar_product(point_of_view, t.get_normal()) > 0.0) {
        if ( true ) {
            os << "\\draw[" << option << "]\n";
            write_face( faceDataPtr, os );
        } // end if
    }     // end for
}

static void myPCG( AMP::LinearAlgebra::Vector::shared_ptr rhs,
                   AMP::LinearAlgebra::Vector::shared_ptr sol,
                   AMP::Operator::Operator::shared_ptr op,
                   AMP::shared_ptr<AMP::Solver::SolverStrategy> pre,
                   size_t maxIters,
                   double relTol,
                   double absTol,
                   bool verbose     = false,
                   std::ostream &os = std::cout )
{
    AMP::LinearAlgebra::Vector::shared_ptr res    = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr dir    = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr ext    = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr oldSol = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr oldRes = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr oldDir = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr matVec = sol->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    op->apply( nullVec, sol, matVec, 1.0, 0.0 );
    oldRes->subtract( rhs, matVec );
    pre->solve( oldRes, ext );
    oldDir->copyVector( ext );
    oldSol->copyVector( sol );
    double initialResNorm = oldRes->L2Norm();
    AMP::AMP_MPI comm     = sol->getComm();
    int rank              = comm.getRank();
    verbose               = verbose && !rank;
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
    AMP::Discretization::DOFManager::shared_ptr dofManager =
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

    // Get the mechanics material model for the contact operator
    AMP::shared_ptr<AMP::Database> model_db = input_db->getDatabase( "MechanicsMaterialModel" );
    AMP::shared_ptr<AMP::Operator::MechanicsModelParameters> mechanicsMaterialModelParams(
        new AMP::Operator::MechanicsModelParameters( model_db ) );
    AMP::shared_ptr<AMP::Operator::MechanicsMaterialModel> masterMechanicsMaterialModel(
        new AMP::Operator::MechanicsMaterialModel( mechanicsMaterialModelParams ) );

    // Build the contact operator
    AMP_INSIST( input_db->keyExists( "ContactOperator" ), "Key ''ContactOperator'' is missing!" );
    AMP::shared_ptr<AMP::Database> contact_db = input_db->getDatabase( "ContactOperator" );
    AMP::shared_ptr<AMP::Operator::ContactOperatorParameters> contactOperatorParams(
        new AMP::Operator::ContactOperatorParameters( contact_db ) );
    contactOperatorParams->d_DOFsPerNode                  = dofsPerNode;
    contactOperatorParams->d_DOFManager                   = dofManager;
    contactOperatorParams->d_GlobalComm                   = globalComm;
    contactOperatorParams->d_Mesh                         = meshAdapter;
    contactOperatorParams->d_MasterMechanicsMaterialModel = masterMechanicsMaterialModel;
    contactOperatorParams->reset(); // got segfault at constructor since d_Mesh was pointing to NULL

    AMP::shared_ptr<AMP::Operator::NodeToGeomType::FaceContactOperator> contactOperator(
        new AMP::Operator::NodeToGeomType::FaceContactOperator( contactOperatorParams ) );

    contactOperator->initialize();

    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> masterBVPOperator;

    // Build the master and slave operators
    AMP::Mesh::MeshID masterMeshID                = contactOperator->getMasterMeshID();
    AMP::Mesh::Mesh::shared_ptr masterMeshAdapter = meshAdapter->Subset( masterMeshID );
    if ( masterMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> masterElementPhysicsModel;
        masterBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                masterMeshAdapter, "MasterBVPOperator", input_db, masterElementPhysicsModel ) );
        columnOperator->append( masterBVPOperator );

        AMP::shared_ptr<AMP::Database> masterSolver_db =
            columnPreconditioner_db->getDatabase( "MasterSolver" );
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> masterSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( masterSolver_db ) );
        masterSolverParams->d_pOperator = masterBVPOperator;
        masterSolverParams->d_comm      = masterMeshAdapter->getComm();
        //    masterSolverParams->d_comm = globalComm;
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> masterSolver(
            new AMP::Solver::PetscKrylovSolver( masterSolverParams ) );
        columnPreconditioner->append( masterSolver );

        std::fstream masterFout;
        masterFout.open( "master_pellet", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 0, masterFout, point_of_view, "blue" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 1, masterFout, point_of_view, "green" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 2, masterFout, point_of_view, "red" );
        drawGeomType::FacesOnBoundaryID(
            masterMeshAdapter, 3, masterFout, point_of_view, "magenta" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 4, masterFout, point_of_view, "black" );
        drawGeomType::FacesOnBoundaryID(
            masterMeshAdapter, 5, masterFout, point_of_view, "orange" );
        drawGeomType::FacesOnBoundaryID( masterMeshAdapter, 6, masterFout, point_of_view, "pink" );
        drawGeomType::FacesOnBoundaryID(
            masterMeshAdapter, 7, masterFout, point_of_view, "violet" );
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 1, masterFout, point_of_view);
        // drawGeomType::FacesOnBoundaryID(masterMeshAdapter, 4, masterFout, point_of_view);
        masterFout.close();
    } // end if

    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> slaveLoadOperator;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> slaveBVPOperator;

    AMP::Mesh::MeshID slaveMeshID                = contactOperator->getSlaveMeshID();
    AMP::Mesh::Mesh::shared_ptr slaveMeshAdapter = meshAdapter->Subset( slaveMeshID );
    if ( slaveMeshAdapter.get() != NULL ) {
        AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> slaveElementPhysicsModel;

        AMP::shared_ptr<AMP::Database> slaveSolver_db =
            columnPreconditioner_db->getDatabase( "SlaveSolver" );
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> slaveSolverParams(
            new AMP::Solver::PetscKrylovSolverParameters( slaveSolver_db ) );

        //    bool useSlaveBVPOperator = input_db->getBool("useSlaveBVPOperator");
        //    if (useSlaveBVPOperator) {
        slaveBVPOperator = AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveBVPOperator", input_db, slaveElementPhysicsModel ) );
        columnOperator->append( slaveBVPOperator );
        slaveSolverParams->d_pOperator = slaveBVPOperator;
        //    } else {
        //      AMP::shared_ptr<AMP::Operator::MechanicsLinearFEOperator>
        //      slaveMechanicsLinearFEOperator =
        //      AMP::dynamic_pointer_cast<
        //          AMP::Operator::MechanicsLinearFEOperator>(AMP::Operator::OperatorBuilder::createOperator(slaveMeshAdapter,
        //                                                                                                   "MechanicsLinearFEOperator",
        //                                                                                                   input_db,
        //                                                                                                   slaveElementPhysicsModel));
        //      columnOperator->append(slaveMechanicsLinearFEOperator);
        //
        //      slaveSolverParams->d_pOperator = slaveMechanicsLinearFEOperator;

        slaveLoadOperator = AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                slaveMeshAdapter, "SlaveLoadOperator", input_db, slaveElementPhysicsModel ) );
        //    AMP::LinearAlgebra::Variable::shared_ptr slaveVar =
        //    slaveMechanicsLinearFEOperator->getOutputVariable();
        AMP::LinearAlgebra::Variable::shared_ptr slaveVar = slaveBVPOperator->getOutputVariable();
        slaveLoadOperator->setVariable( slaveVar );
        //    } // end if

        //    slaveSolverParams->d_comm = globalComm;
        slaveSolverParams->d_comm = slaveMeshAdapter->getComm();
        AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> slaveSolver(
            new AMP::Solver::PetscKrylovSolver( slaveSolverParams ) );
        columnPreconditioner->append( slaveSolver );

        std::fstream slaveFout;
        slaveFout.open( "slave_pellet", std::fstream::out );
        double point_of_view[3] = { 1.0, 1.0, 1.0 };
        drawGeomType::FacesOnBoundaryID(
            slaveMeshAdapter, 0, slaveFout, point_of_view, "dashed,red" );
        // drawGeomType::FacesOnBoundaryID(slaveMeshAdapter, 1, slaveFout, point_of_view, "dashed");
        // drawGeomType::FacesOnBoundaryID(slaveMeshAdapter, 4, slaveFout, point_of_view, "dashed");
        // drawVerticesOnBoundaryID(slaveMeshAdapter, 2, slaveFout, point_of_view, "red");
        slaveFout.close();
    } // end if

    AMP::shared_ptr<AMP::Database> contactPreconditioner_db =
        columnPreconditioner_db->getDatabase( "ContactPreconditioner" );
    AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolverParameters>
        contactPreconditionerParams(
            new AMP::Solver::ConstraintsEliminationSolverParameters( contactPreconditioner_db ) );
    contactPreconditionerParams->d_pOperator = contactOperator;
    AMP::shared_ptr<AMP::Solver::ConstraintsEliminationSolver> contactPreconditioner(
        new AMP::Solver::ConstraintsEliminationSolver( contactPreconditionerParams ) );
    columnPreconditioner->append( contactPreconditioner );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Variable::shared_ptr columnVar = columnOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr columnSolVec =
        createVector( dofManager, columnVar, split );
    AMP::LinearAlgebra::Vector::shared_ptr columnRhsVec =
        createVector( dofManager, columnVar, split );
    columnSolVec->zero();
    columnRhsVec->zero();

    bool skipDisplaceMesh = true;
    contactOperator->updateActiveSet( nullVec, skipDisplaceMesh );

#ifdef USE_EXT_SILO
    {
        siloWriter->registerVector(
            columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
        char outFileName[256];
        sprintf( outFileName, "TOTO_%d", 0 );
        siloWriter->writeFile( outFileName, 0 );
    }
#endif

    size_t const maxActiveSetIterations =
        input_db->getIntegerWithDefault( "maxActiveSetIterations", 5 );
    for ( size_t activeSetIteration = 0; activeSetIteration < maxActiveSetIterations;
          ++activeSetIteration ) {
        if ( !rank ) {
            std::cout << "ACTIVE SET ITERATION #" << activeSetIteration + 1 << "\n";
        }

        columnSolVec->zero();
        columnRhsVec->zero();

        //  if (slaveLoadOperator.get() != NULL) {
        //    slaveLoadOperator->apply(nullVec, nullVec, columnRhsVec, 1.0, 0.0);
        //  } // end if

        // apply dirichlet rhs correction
        if ( masterBVPOperator.get() != NULL ) {
            masterBVPOperator->modifyRHSvector( columnRhsVec );
        } // end if
        if ( slaveBVPOperator.get() != NULL ) {
            slaveBVPOperator->modifyRHSvector( columnRhsVec );
        } // end if

        // get d
        contactOperator->addShiftToSlave( columnSolVec );

        // compute - Kd
        AMP::LinearAlgebra::Vector::shared_ptr rhsCorrectionVec =
            createVector( dofManager, columnVar, split );
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

        bool usePetscKrylovSolver = input_db->getBool( "usePetscKrylovSolver" );
        if ( usePetscKrylovSolver ) {
            // Build a matrix shell operator to use the column operator with the petsc krylov
            // solvers
            AMP::shared_ptr<AMP::Database> matrixShellDatabase =
                input_db->getDatabase( "MatrixShellOperator" );
            AMP::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(
                new AMP::Operator::OperatorParameters( matrixShellDatabase ) );
            AMP::shared_ptr<AMP::Operator::PetscMatrixShellOperator> matrixShellOperator(
                new AMP::Operator::PetscMatrixShellOperator( matrixShellParams ) );

            int numMasterLocalNodes = 0;
            int numSlaveLocalNodes  = 0;
            if ( masterMeshAdapter.get() != NULL ) {
                numMasterLocalNodes =
                    masterMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
            }
            if ( slaveMeshAdapter.get() != NULL ) {
                numSlaveLocalNodes =
                    slaveMeshAdapter->numLocalElements( AMP::Mesh::GeomType::Vertex );
            }
            int matLocalSize = dofsPerNode * ( numMasterLocalNodes + numSlaveLocalNodes );
            AMP_ASSERT( matLocalSize == static_cast<int>( dofManager->numLocalDOF() ) );
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

            linearSolver->solve( columnRhsVec, columnSolVec );
        } else {
            size_t myPCGmaxIters = input_db->getInteger( "myPCGmaxIters" );
            double myPCGrelTol   = input_db->getDouble( "myPCGrelTol" );
            double myPCGabsTol   = input_db->getDouble( "myPCGabsTol" );
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

#ifdef USE_EXT_SILO
        meshAdapter->displaceMesh( columnSolVec );
        siloWriter->registerVector(
            columnSolVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
        char outFileName[256];
        sprintf( outFileName, "TOTO_%d", 0 );
        siloWriter->writeFile( outFileName, activeSetIteration + 1 );
        columnSolVec->scale( -1.0 );
        meshAdapter->displaceMesh( columnSolVec );
        columnSolVec->scale( -1.0 );
#endif

        //  meshAdapter->displaceMesh(columnSolVec);
        size_t nChangesInActiveSet = contactOperator->updateActiveSet( columnSolVec );
        if ( !rank ) {
            std::cout << nChangesInActiveSet << " CHANGES IN ACTIVE SET\n";
        }

        if ( nChangesInActiveSet == 0 ) {
            break;
        }
        AMP_ASSERT( activeSetIteration != maxActiveSetIterations - 1 );
    } // end for
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


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.push_back( "testNodeToSegmentConstraintsOperator-cube" );
    //  exeNames.push_back("testNodeToSegmentConstraintsOperator-cylinder");
    //  exeNames.push_back("testNodeToSegmentConstraintsOperator-pellet");

    for ( size_t i = 0; i < exeNames.size(); ++i ) {
        myTest( &ut, exeNames[i] );
    } // end for

    ut.report();
    int num_failed = ut.NumFailGlobal();

    AMP::AMPManager::shutdown();
    return num_failed;
}
