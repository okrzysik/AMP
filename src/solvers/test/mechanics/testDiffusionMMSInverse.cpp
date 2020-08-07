#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/operators/BVPOperatorParameters.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NonlinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrectionParameters.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionNonlinearElement.h"
#include "AMP/operators/diffusion/DiffusionNonlinearFEOperator.h"
#include "AMP/operators/libmesh/MassDensityModel.h"
#include "AMP/operators/libmesh/MassLinearFEOperator.h"
#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/ManufacturedSolution.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include "../../../operators/test/applyTests.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <string>


static void inverseTest1( AMP::UnitTest *ut, const std::string &exeName )
{
    // Tests diffusion Dirchlet BVP operator for temperature

    // Initialization
    std::string input_file = exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Get the Mesh database and create the mesh parameters
    std::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> params( new AMP::Mesh::MeshParameters( database ) );
    params->setComm( globalComm );

    // Create the meshes from the input database
    AMP::Mesh::Mesh::shared_ptr meshAdapter = AMP::Mesh::Mesh::buildMesh( params );

    // Create nonlinear diffusion BVP operator and access volume nonlinear Diffusion operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> nonlinearPhysicsModel;
    std::shared_ptr<AMP::Operator::Operator> nlinBVPOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "FickNonlinearBVPOperator", input_db, nonlinearPhysicsModel );
    std::shared_ptr<AMP::Operator::NonlinearBVPOperator> nlinBVPOp =
        std::dynamic_pointer_cast<AMP::Operator::NonlinearBVPOperator>( nlinBVPOperator );
    std::shared_ptr<AMP::Operator::DiffusionNonlinearFEOperator> nlinOp =
        std::dynamic_pointer_cast<AMP::Operator::DiffusionNonlinearFEOperator>(
            nlinBVPOp->getVolumeOperator() );

    // Acquire Dirichlet boundary operator and parameters
    std::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletOp =
        std::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            nlinBVPOp->getBoundaryOperator() );
    // std::shared_ptr<AMP::Database> dirichlet_db =
    // input_db->getDatabase("FickDirichletVectorCorrection");
    // AMP_INSIST(dirichlet_db->getScalar<int>("valuesType")==2,
    // "DirichletVectorCorrection::valuesType must be 2");
    // dirichlet_db->putScalar("isAttachedToVolumeOperator",false); // OperatorBuilder forces this
    // to be true for some reason.
    // std::shared_ptr<AMP::Operator::DirichletVectorCorrectionParameters> dirichletParams(
    //      new AMP::Operator::DirichletVectorCorrectionParameters(dirichlet_db));
    // std::shared_ptr<AMP::Operator::OperatorParameters> dirichletOpParams =
    // std::dynamic_pointer_cast<AMP::Operator::OperatorParameters>(
    //      dirichletParams);
    // dirichletOp->reset(dirichletOpParams);

    // Create linear diffusion BVP operator with bc's
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> linearPhysicsModel;
    std::shared_ptr<AMP::Operator::Operator> linBVPOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "FickLinearBVPOperator", input_db, linearPhysicsModel );
    std::shared_ptr<AMP::Operator::LinearBVPOperator> linBVPOp =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>( linBVPOperator );

    // Get source mass operator
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> sourcePhysicsModel;
    std::shared_ptr<AMP::Operator::Operator> sourceOperator =
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "ManufacturedSourceOperator", input_db, sourcePhysicsModel );
    std::shared_ptr<AMP::Operator::MassLinearFEOperator> sourceOp =
        std::dynamic_pointer_cast<AMP::Operator::MassLinearFEOperator>( sourceOperator );

    std::shared_ptr<AMP::Operator::MassDensityModel> densityModel = sourceOp->getDensityModel();
    std::shared_ptr<AMP::ManufacturedSolution> mfgSolution =
        densityModel->getManufacturedSolution();

    // Set up input and output vectors
    AMP::LinearAlgebra::Variable::shared_ptr solVar = nlinOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr rhsVar = nlinOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr resVar = nlinOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr bndVar = nlinOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr inpVar = sourceOp->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr srcVar = sourceOp->getOutputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr workVar( new AMP::LinearAlgebra::Variable( "work" ) );

    AMP::Discretization::DOFManager::shared_ptr DOF = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 1, true );

    AMP::LinearAlgebra::Vector::shared_ptr solVec = AMP::LinearAlgebra::createVector( DOF, solVar );
    AMP::LinearAlgebra::Vector::shared_ptr rhsVec = AMP::LinearAlgebra::createVector( DOF, rhsVar );
    AMP::LinearAlgebra::Vector::shared_ptr resVec = AMP::LinearAlgebra::createVector( DOF, resVar );
    AMP::LinearAlgebra::Vector::shared_ptr bndVec = AMP::LinearAlgebra::createVector( DOF, bndVar );
    AMP::LinearAlgebra::Vector::shared_ptr inpVec = AMP::LinearAlgebra::createVector( DOF, inpVar );
    AMP::LinearAlgebra::Vector::shared_ptr srcVec = AMP::LinearAlgebra::createVector( DOF, srcVar );
    AMP::LinearAlgebra::Vector::shared_ptr workVec =
        AMP::LinearAlgebra::createVector( DOF, workVar );

    srcVec->setToScalar( 0., srcVec );

    // Fill in manufactured solution in mesh interior
    const double Pi                  = 3.1415926535898;
    AMP::Mesh::MeshIterator iterator = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
    std::string mfgName              = mfgSolution->get_name();
    bool isCylindrical               = mfgName.find( "Cylindrical" ) < mfgName.size();
    for ( ; iterator != iterator.end(); ++iterator ) {
        double x, y, z;
        std::valarray<double> poly( 10 );
        auto coord = iterator->coord();
        x          = coord[0];
        y          = coord[1];
        z          = coord[2];
        if ( isCylindrical ) {
            double th = 0.;
            double r  = sqrt( x * x + y * y );
            if ( r > 0 ) {
                th = acos( x / r );
                if ( y < 0. )
                    th = 2 * Pi - th;
            }
            mfgSolution->evaluate( poly, r, th, z );
        } else {
            mfgSolution->evaluate( poly, x, y, z );
        }
        std::vector<size_t> gid;
        DOF->getDOFs( iterator->globalID(), gid );
        std::vector<double> srcVal( 1 ), dumT( 1 ), dumU( 1 ), dumB( 1 );
        std::vector<libMesh::Point> point( 1, libMesh::Point( x, y, z ) );
        densityModel->getDensityManufactured( srcVal, dumT, dumU, dumB, point );
        inpVec->setValueByGlobalID( gid[0], srcVal[0] );
    }

    // Fill in manufactured solution on mesh boundary
    for ( int j = 0; j <= 8; j++ ) {
        AMP::Mesh::MeshIterator beg_bnd =
            meshAdapter->getBoundaryIDIterator( AMP::Mesh::GeomType::Vertex, j, 0 );
        AMP::Mesh::MeshIterator end_bnd = beg_bnd.end();
        AMP::Mesh::MeshIterator iter;
        for ( iter = beg_bnd; iter != end_bnd; ++iter ) {
            std::valarray<double> poly( 10 );
            double x, y, z;
            auto coord = iterator->coord();
            x          = coord[0];
            y          = coord[1];
            z          = coord[2];
            if ( isCylindrical ) {
                double th = 0.;
                double r  = sqrt( x * x + y * y );
                if ( r > 0 ) {
                    th = acos( x / r );
                    if ( y < 0. )
                        th = 2 * Pi - th;
                }
                mfgSolution->evaluate( poly, r, th, z );
            } else {
                mfgSolution->evaluate( poly, x, y, z );
            }
            std::vector<size_t> gid;
            DOF->getDOFs( iterator->globalID(), gid );
            bndVec->setValueByGlobalID( gid[0], poly[0] );
        }
    }

    // Set boundary values for manufactured solution for sinusoid, gaussian, etc. (non constant BC)
    dirichletOp->setVariable( bndVar );
    dirichletOp->setDirichletValues( bndVec );

    // Evaluate manufactured solution as an FE source
    sourceOp->apply( inpVec, rhsVec );

    // Reset solution vector to initial value and print out norm
    solVec->setToScalar( 0.1, solVec );

    // Set up initial guess
    nlinBVPOp->modifyInitialSolutionVector( solVec );

    // Set up solver
    std::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    std::shared_ptr<AMP::Database> linearSolver_db =
        nonlinearSolver_db->getDatabase( "LinearSolver" );
    std::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = nlinBVPOp;
    nonlinearSolverParams->d_pInitialGuess = solVec;
    std::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    // Set up preconditioner
    std::shared_ptr<AMP::Database> preconditioner_db =
        linearSolver_db->getDatabase( "Preconditioner" );
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> preconditionerParams(
        new AMP::Solver::SolverStrategyParameters( preconditioner_db ) );
    preconditionerParams->d_pOperator = linBVPOp;
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> preconditioner(
        new AMP::Solver::TrilinosMLSolver( preconditionerParams ) );

    // Register the preconditioner with the Jacobian free Krylov solver
    std::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();
    linearSolver->setPreconditioner( preconditioner );

    // Get initial residual
    nlinBVPOp->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm( resVec );
    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    // Run solver
    nonlinearSolver->setZeroInitialGuess( false );
    nonlinearSolver->solve( rhsVec, solVec );

    // Get final residual
    nlinBVPOp->residual( rhsVec, solVec, resVec );
    double finalResidualNorm = resVec->L2Norm( resVec );
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    // Final communication
    solVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    resVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    // Output Mathematica form (requires serial execution)
    for ( int i = 0; i < globalComm.getSize(); i++ ) {
        if ( globalComm.getRank() == i ) {
            std::string filename          = "data_" + exeName;
            int rank                      = globalComm.getRank();
            int nranks                    = globalComm.getSize();
            std::ios_base::openmode omode = std::ios_base::out;
            if ( rank > 0 )
                omode |= std::ios_base::app;
            std::ofstream file( filename.c_str(), omode );
            if ( rank == 0 ) {
                file << "(* x y z solution solution fe-source fe-operator error *)" << std::endl;
                file << "results={" << std::endl;
            }

            iterator        = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
            size_t numNodes = 0;
            for ( ; iterator != iterator.end(); ++iterator )
                numNodes++;

            iterator     = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, 0 );
            size_t iNode = 0;
            double l2err = 0.;
            for ( ; iterator != iterator.end(); ++iterator ) {
                double x, y, z;
                auto coord = iterator->coord();
                x          = coord[0];
                y          = coord[1];
                z          = coord[2];
                std::vector<size_t> gid;
                DOF->getDOFs( iterator->globalID(), gid );
                double val, res, sol, src, err;
                res = resVec->getValueByGlobalID( gid[0] );
                sol = solVec->getValueByGlobalID( gid[0] );
                src = srcVec->getValueByGlobalID( gid[0] );
                err = res / ( src + .5 * res + std::numeric_limits<double>::epsilon() );
                std::valarray<double> poly( 10 );
                mfgSolution->evaluate( poly, x, y, z );
                val = poly[0];
                workVec->setValueByGlobalID( gid[0], err );

                file << "{" << x << "," << y << "," << z << "," << val << "," << sol << "," << src
                     << "," << res + src << "," << err << "}";
                if ( iNode < numNodes - 1 )
                    file << "," << std::endl;

                l2err += ( res * res );
                iNode++;
            }

            if ( rank == nranks - 1 ) {
                file << "};" << std::endl;
                file << "nodes = " << numNodes << "; l2err = " << l2err << ";" << std::endl;
            }

            file.close();
        }
        globalComm.barrier();
    }

    // Plot the results
    if ( globalComm.getSize() == 1 ) {
#ifdef USE_EXT_SILO
        AMP::Utilities::Writer::shared_ptr siloWriter =
            AMP::Utilities::Writer::buildWriter( "Silo" );
        siloWriter->registerVector(
            workVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "RelativeError" );
        siloWriter->registerVector( solVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Solution" );
        siloWriter->registerVector( srcVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Source" );
        siloWriter->registerVector( resVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );
        siloWriter->writeFile( exeName, 0 );
#endif
    }

    ut->passes( exeName );
    std::cout.flush();
}

int testDiffusionMMSInverse( int argc, char *argv[] )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    // Check to see if an input file was requested on the command line
    std::vector<std::string> files;
    std::vector<std::string> arguments( argv + 1, argv + argc );
    if ( argc > 1 ) {
        // Populate array with argv - easier with which to work
        for ( size_t i = 0; i < arguments.size(); ++i ) {
            if ( arguments[i][0] == '-' )
                i++; // Move past the next argument - not a filename
            else
                files.push_back( arguments[i] ); // Store this as a file
        }
    } else {
        std::cout
            << "No input files are currently hardcoded. Files must be given as an argument.\n";
        return 1;
    }

    // Run the tests
    for ( auto &file : files ) {
        inverseTest1( &ut, file );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
