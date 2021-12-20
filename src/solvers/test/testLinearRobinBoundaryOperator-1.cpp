#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/MeshParameters.h"

#include "AMP/IO/PIO.h"
#include "AMP/IO/Writer.h"
#include "AMP/operators/ElementOperationFactory.h"
#include "AMP/operators/ElementPhysicsModelFactory.h"
#include "AMP/operators/LinearBVPOperator.h"
#include "AMP/operators/OperatorBuilder.h"
#include "AMP/operators/boundary/ColumnBoundaryOperator.h"
#include "AMP/operators/boundary/DirichletMatrixCorrection.h"
#include "AMP/operators/boundary/DirichletVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/NeumannVectorCorrection.h"
#include "AMP/operators/boundary/libmesh/RobinMatrixCorrection.h"
#include "AMP/operators/diffusion/DiffusionLinearElement.h"
#include "AMP/operators/diffusion/DiffusionLinearFEOperator.h"
#include "AMP/operators/diffusion/DiffusionTransportModel.h"
#include "AMP/operators/libmesh/VolumeIntegralOperator.h"
#include "AMP/solvers/trilinos/ml/TrilinosMLSolver.h"
#include "AMP/utils/AMPManager.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Database.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"

#include <memory>
#include <string>


constexpr double pi = 3.14159265;
static inline double fun_T0( double x, double y, double z, double sg )
{
    return ( sg * cos( 0.1 * pi * x ) * cos( 0.1 * pi * y ) * cos( 0.1 * pi * z ) );
}
static inline double fun_dTdx( double x, double y, double z, double sg )
{
    return ( sg * -0.1 * pi * sin( 0.1 * pi * x ) * cos( 0.1 * pi * y ) * cos( 0.1 * pi * z ) );
}
static inline double fun_dTdy( double x, double y, double z, double sg )
{
    return ( sg * -0.1 * pi * cos( 0.1 * pi * x ) * sin( 0.1 * pi * y ) * cos( 0.1 * pi * z ) );
}
static inline double fun_dTdz( double x, double y, double z, double sg )
{
    return ( sg * -0.1 * pi * cos( 0.1 * pi * x ) * cos( 0.1 * pi * y ) * sin( 0.1 * pi * z ) );
}
static inline double fun_rhs( double x, double y, double z, double sg )
{
    return ( sg * -0.03 * pi * pi * cos( 0.1 * pi * x ) * cos( 0.1 * pi * y ) *
             cos( 0.1 * pi * z ) );
}


void linearRobinTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Input and output file names
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;

    // Fill the database from the input file.
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );

    // Print from all cores into the output files
    AMP::logAllNodes( log_file );

    //--------------------------------------------------
    //   Create the Mesh.
    //--------------------------------------------------
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    auto mesh_db   = input_db->getDatabase( "Mesh" );
    auto mgrParams = std::make_shared<AMP::Mesh::MeshParameters>( mesh_db );
    mgrParams->setComm( AMP::AMP_MPI( AMP_COMM_WORLD ) );
    auto meshAdapter = AMP::Mesh::Mesh::buildMesh( mgrParams );

    //--------------------------------------------------
    //   Create DOF Managers.
    //--------------------------------------------------
    int DOFsPerElement = 8;
    int DOFsPerNode    = 1;
    int ghostWidth     = 1;
    bool split         = true;
    auto nodalDofMap   = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Vertex, ghostWidth, DOFsPerNode, split );
    auto gaussPointDofMap = AMP::Discretization::simpleDOFManager::create(
        meshAdapter, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerElement, split );

    // Create a shared pointer to a Variable - Power - Output because it will be used in the
    // "residual" location of apply.
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    //------------------------------------------
    //   CREATE THE THERMAL BVP OPERATOR  //
    //------------------------------------------
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    auto diffusionOperator = std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );

    auto TemperatureInKelvinVec = AMP::LinearAlgebra::createVector(
        nodalDofMap, diffusionOperator->getOutputVariable(), split );
    auto RightHandSideVec = TemperatureInKelvinVec->cloneVector();
    auto ResidualVec      = TemperatureInKelvinVec->cloneVector();
    auto variableFluxVec  = TemperatureInKelvinVec->cloneVector();

    RightHandSideVec->zero();
    variableFluxVec->zero();

    //------------------------------------------

    auto boundaryOp = diffusionOperator->getBoundaryOperator();

    auto robinBoundaryOp =
        std::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>( boundaryOp )
            ->getBoundaryOperator( 0 );

    auto robinboundaryDatabase = input_db->getDatabase( "RobinMatrixCorrection" );

    robinboundaryDatabase->putScalar( "constant_flux", false );
    robinboundaryDatabase->putScalar( "skip_matrix_correction", true );
    auto correctionParameters =
        std::make_shared<AMP::Operator::RobinMatrixCorrectionParameters>( robinboundaryDatabase );
    //------------------------------------------


    //------------------------------------------
    // check the solution
    int zeroGhostWidth = 0;
    auto node          = meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    auto end_node      = node.end();

    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double val, rhs;

        rhs = fun_rhs( px, py, pz, -1.0 );
        RightHandSideVec->setValuesByGlobalID( 1, &gid[0], &rhs );

        if ( fabs( pz - 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdz( px, py, pz, 1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        } else if ( fabs( pz + 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdz( px, py, pz, -1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        } else if ( fabs( px - 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdx( px, py, pz, 1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        } else if ( fabs( px + 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdx( px, py, pz, -1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        } else if ( fabs( py - 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdy( px, py, pz, 1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        } else if ( fabs( py + 1.0 ) <= 1.0e-12 ) {
            val = fun_dTdy( px, py, pz, -1.0 );
            val = val + fun_T0( px, py, pz, 1.0 );
            variableFluxVec->setValuesByGlobalID( 1, &gid[0], &val );
        }
    } // end for node
    RightHandSideVec->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    correctionParameters->d_variableFlux = variableFluxVec;
    robinBoundaryOp->reset( correctionParameters );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    auto sourceOperator = std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
        AMP::Operator::OperatorBuilder::createOperator(
            meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    auto SourceVar = sourceOperator->getOutputVariable();
    auto SourceVec = AMP::LinearAlgebra::createVector( nodalDofMap, SourceVar, split );
    SourceVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( RightHandSideVec, SourceVec );

    //------------------------------------------
    //   Add the boundary conditions corrections //
    //------------------------------------------

    std::cout << "RHS Norm before BC Correction " << SourceVec->L2Norm() << std::endl;

    diffusionOperator->modifyRHSvector( SourceVec );

    auto rhsNorm = static_cast<double>( SourceVec->L2Norm() );
    std::cout << "RHS Norm after BC Correction " << rhsNorm << std::endl;

    //------------------------------------------

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    auto mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters fo the class with the info on the database.
    auto mlSolverParams = std::make_shared<AMP::Solver::SolverStrategyParameters>( mlSolver_db );

    // Define the operature to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = static_cast<double>( TemperatureInKelvinVec->L2Norm() );
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    // Create the ML Solver
    auto mlSolver = std::make_shared<AMP::Solver::TrilinosMLSolver>( mlSolverParams );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    // Solve the prblem.
    mlSolver->apply( SourceVec, TemperatureInKelvinVec );

    // Compute the residual
    diffusionOperator->residual( SourceVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = static_cast<double>( ResidualVec->L2Norm() );
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    node          = node.begin();
    auto diffVec  = TemperatureInKelvinVec->cloneVector();
    auto exactVec = TemperatureInKelvinVec->cloneVector();

    diffVec->zero();
    exactVec->zero();

    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double exact;
        exact = fun_T0( px, py, pz, 1.0 );
        exactVec->setValuesByGlobalID( 1, &gid[0], &exact );
    }

    diffVec->subtract( *exactVec, *TemperatureInKelvinVec );

    double exactNorm = static_cast<double>( exactVec->L1Norm() );
    std::cout << "L2norm of exactVec " << exactNorm << std::endl;

    double solutionNorm = static_cast<double>( TemperatureInKelvinVec->L1Norm() );
    std::cout << "L2norm of solutionVec " << solutionNorm << std::endl;

    double errorNorm = static_cast<double>( diffVec->L1Norm() );
    std::cout << "L1norm of DiffVec " << errorNorm << std::endl;

    if ( errorNorm > 1.0 ) {
        ut->failure( "linear robin boundary operator verification test-1." );
    } else {
        ut->passes( "linear robin boundary operator verification test-1." );
    }

    // Plot the results
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    auto siloWriter = AMP::IO::Writer::buildWriter( "Silo" );
    siloWriter->registerVector(
        TemperatureInKelvinVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "TemperatureInKelvin" );
    siloWriter->registerVector( exactVec, meshAdapter, AMP::Mesh::GeomType::Vertex, "Exact" );
    siloWriter->writeFile( input_file, 0 );
#endif

    ut->passes( exeName );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;
    exeNames.emplace_back( "testLinearRobinBoundaryOperator-1" );
    // exeNames.push_back("testLinearRobinBoundaryOperator-2");

    for ( auto &exeName : exeNames ) {
        linearRobinTest( &ut, exeName );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
