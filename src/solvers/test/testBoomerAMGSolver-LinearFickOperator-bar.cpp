#include "AMP/ampmesh/Mesh.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/NeutronicsRhs.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/hypre/BoomerAMGSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <fstream>
#include <limits>
#include <string>

#define ITFAILS ut.failure( __LINE__ );
#define UNIT_TEST( a ) \
    if ( !( a ) )      \
        ut.failure( __LINE__ );

void linearFickTest( AMP::UnitTest *ut )
{
    // Input and output file names
    std::string exeName( "testBoomerAMGSolver-LinearFickOperator-bar" );
    std::string input_file  = "input_" + exeName;
    std::string log_file    = "output_" + exeName;
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::PIO::logAllNodes( log_file );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    std::shared_ptr<AMP::Database> mesh_db = input_db->getDatabase( "Mesh" );
    std::shared_ptr<AMP::Mesh::MeshParameters> mgrParams(
        new AMP::Mesh::MeshParameters( mesh_db ) );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    std::shared_ptr<AMP::Mesh::Mesh> meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );
    //--------------------------------------------------

    //--------------------------------------------------
    // Create a DOF manager for a nodal vector
    //--------------------------------------------------
    int DOFsPerNode     = 1;
    int nodalGhostWidth = 1;
    bool split          = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, nodalGhostWidth, DOFsPerNode, split );
    //--------------------------------------------------

    ////////////////////////////////////
    //   CREATE THE DIFFUSION OPERATOR  //
    ////////////////////////////////////

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );

    AMP::LinearAlgebra::Vector::shared_ptr SolutionVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getInputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, diffusionOperator->getOutputVariable() );

    RightHandSideVec->setToScalar( 0. );

    std::shared_ptr<AMP::Operator::BoundaryOperator> boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    boundaryOp->addRHScorrection( RightHandSideVec );
    boundaryOp->setRHScorrection( RightHandSideVec );

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    std::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters for the class with the info on the database.
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
        new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

    // Define the operator to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    //////////////////////////
    //   FIND THE SOLUTION  //
    //////////////////////////

    // Set initial guess
    SolutionVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = SolutionVec->L2Norm();
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    double rhsNorm = RightHandSideVec->L2Norm();
    std::cout << "RHS Norm: " << rhsNorm << std::endl;

    // Create the ML Solver
    auto mlSolver = std::make_shared<AMP::Solver::BoomerAMGSolver>( mlSolverParams );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    // Solve the problem.
    mlSolver->solve( RightHandSideVec, SolutionVec );

    // Compute the residual
    diffusionOperator->residual( RightHandSideVec, SolutionVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm();
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > 10.0 ) {
        ut->failure( "BoomerAMGSolver unsuccessfully solves a linear fick problem." );
    } else {
        ut->passes( "BoomerAMGSolver successfully solves a linear fick problem." );
    }

    ///////////////////////////
    //   CHECK THE SOLUTION  //
    ///////////////////////////

    int zeroGhostWidth = 0;
    AMP::Mesh::MeshIterator iterator =
        meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );

    // The analytical solution is:  T = a + b*z + c*z*z
    //   c = -power/2
    //   b = -10*power
    //   a = 300 + 150*power

    double power = 1.;
    double c     = -power / 2.;
    double b     = -10. * power;
    double a     = 300. + 150. * power;
    bool passes  = true;

    if ( false ) {
        double cal, zee, sol, err;
        // Serialize the code
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
                    file << "(* x y z analytic calculated relative-error *)" << std::endl;
                    file << "formula=" << a << " + " << b << "*z + " << c << "*z^2;" << std::endl;
                    file << "results={" << std::endl;
                }
                file.precision( 14 );

                iterator        = iterator.begin();
                size_t numNodes = 0, iNode = 0;
                for ( ; iterator != iterator.end(); ++iterator )
                    numNodes++;

                iterator = iterator.begin();
                for ( ; iterator != iterator.end(); ++iterator ) {
                    std::vector<size_t> gid;
                    nodalDofMap->getDOFs( iterator->globalID(), gid );
                    cal = SolutionVec->getValueByGlobalID( gid[0] );
                    zee = ( iterator->coord() )[2];
                    sol = a + b * zee + c * zee * zee;
                    err = fabs( cal - sol ) * 2. /
                          ( cal + sol + std::numeric_limits<double>::epsilon() );
                    double x, y, z;
                    x = ( iterator->coord() )[0];
                    y = ( iterator->coord() )[1];
                    z = ( iterator->coord() )[2];
                    file << "{" << x << "," << y << "," << z << "," << sol << "," << cal << ","
                         << err << "}";
                    if ( iNode < numNodes - 1 )
                        file << "," << std::endl;
                    if ( fabs( cal - sol ) > cal * 1e-3 ) {
                        passes = 0;
                        ut->failure( "Error" );
                    }
                    iNode++;
                }

                if ( rank == nranks - 1 ) {
                    file << "};" << std::endl;
                }
                file.close();
            }
        }
        if ( passes )
            ut->passes( "The linear fick solve is verified." );
    }

// Plot the results
#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
    siloWriter->registerMesh( meshAdapter );

    siloWriter->registerVector(
        SolutionVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Concentration" );
    siloWriter->registerVector( ResidualVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Residual" );

    siloWriter->writeFile( input_file, 0 );
#endif

    input_db.reset();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    linearFickTest( &ut );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
