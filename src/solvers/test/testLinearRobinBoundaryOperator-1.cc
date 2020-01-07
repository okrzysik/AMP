#include "AMP/discretization/DOF_Manager.h"
#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/materials/Material.h"
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
#include "AMP/utils/PIO.h"
#include "AMP/utils/UnitTest.h"
#include "AMP/utils/Utilities.h"
#include "AMP/utils/Writer.h"
#include "AMP/vectors/Variable.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/VectorBuilder.h"
#include <memory>

#include <string>


#define __PI__ 3.14159265
#define __INIT_T0__( x, y, z, sg ) \
    ( sg * cos( 0.1 * __PI__ * x ) * cos( 0.1 * __PI__ * y ) * cos( 0.1 * __PI__ * z ) )
#define __INIT_dTdx__( x, y, z, sg )                                           \
    ( sg * -0.1 * __PI__ * sin( 0.1 * __PI__ * x ) * cos( 0.1 * __PI__ * y ) * \
      cos( 0.1 * __PI__ * z ) )
#define __INIT_dTdy__( x, y, z, sg )                                           \
    ( sg * -0.1 * __PI__ * cos( 0.1 * __PI__ * x ) * sin( 0.1 * __PI__ * y ) * \
      cos( 0.1 * __PI__ * z ) )
#define __INIT_dTdz__( x, y, z, sg )                                           \
    ( sg * -0.1 * __PI__ * cos( 0.1 * __PI__ * x ) * cos( 0.1 * __PI__ * y ) * \
      sin( 0.1 * __PI__ * z ) )
#define __INIT_rhs__( x, y, z, sg )                                                      \
    ( sg * -0.03 * __PI__ * __PI__ * cos( 0.1 * __PI__ * x ) * cos( 0.1 * __PI__ * y ) * \
      cos( 0.1 * __PI__ * z ) )


void linearRobinTest( AMP::UnitTest *ut, const std::string &exeName )
{
    // Input and output file names
    //  #include <string>
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    ////////////////////////////////////
    //    INITIALIZE THE PROBLEM      //
    ////////////////////////////////////

    // Create the map to get an available material from a string.
    //  #include "AMP/materials/Material.h"

    // Construct a smart pointer to a new database.
    //  #include <memory>
    //  #include "AMP/utils/Database.h"


    // Fill the database from the input file.
    //  #include "AMP/utils/Database.h"
    auto input_db = AMP::Database::parseInputFile( input_file );
    input_db->print( AMP::plog );


    // Print from all cores into the output files
    //   #include "AMP/utils/PIO.h"
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
    //   Create DOF Managers.
    //--------------------------------------------------
    int DOFsPerElement = 8;
    int DOFsPerNode    = 1;
    int ghostWidth     = 1;
    bool split         = true;
    AMP::Discretization::DOFManager::shared_ptr nodalDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, ghostWidth, DOFsPerNode, split );
    AMP::Discretization::DOFManager::shared_ptr gaussPointDofMap =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Volume, ghostWidth, DOFsPerElement, split );

    // Create a shared pointer to a Variable - Power - Output because it will be used in the
    // "residual" location of
    // apply.
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    //------------------------------------------
    //   CREATE THE THERMAL BVP OPERATOR  //
    //------------------------------------------
    std::shared_ptr<AMP::Operator::ElementPhysicsModel> transportModel;
    std::shared_ptr<AMP::Operator::LinearBVPOperator> diffusionOperator =
        std::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "DiffusionBVPOperator", input_db, transportModel ) );

    AMP::LinearAlgebra::Vector::shared_ptr TemperatureInKelvinVec =
        AMP::LinearAlgebra::createVector(
            nodalDofMap, diffusionOperator->getOutputVariable(), split );
    AMP::LinearAlgebra::Vector::shared_ptr RightHandSideVec = TemperatureInKelvinVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr ResidualVec      = TemperatureInKelvinVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr variableFluxVec  = TemperatureInKelvinVec->cloneVector();

    RightHandSideVec->zero();
    variableFluxVec->zero();

    //------------------------------------------

    AMP::Operator::Operator::shared_ptr boundaryOp;
    boundaryOp = diffusionOperator->getBoundaryOperator();

    AMP::Operator::Operator::shared_ptr robinBoundaryOp;
    robinBoundaryOp =
        ( std::dynamic_pointer_cast<AMP::Operator::ColumnBoundaryOperator>( boundaryOp ) )
            ->getBoundaryOperator( 0 );

    std::shared_ptr<AMP::Database> robinboundaryDatabase = std::dynamic_pointer_cast<AMP::Database>(
        input_db->getDatabase( "RobinMatrixCorrection" ) );

    robinboundaryDatabase->putScalar( "constant_flux", false );
    robinboundaryDatabase->putScalar( "skip_matrix_correction", true );
    std::shared_ptr<AMP::Operator::RobinMatrixCorrectionParameters> correctionParameters(
        new AMP::Operator::RobinMatrixCorrectionParameters( robinboundaryDatabase ) );
    //------------------------------------------


    //------------------------------------------
    // check the solution
    int zeroGhostWidth = 0;
    AMP::Mesh::MeshIterator node =
        meshAdapter->getIterator( AMP::Mesh::GeomType::Vertex, zeroGhostWidth );
    AMP::Mesh::MeshIterator end_node = node.end();

    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double val, rhs;

        rhs = __INIT_rhs__( px, py, pz, -1.0 );
        RightHandSideVec->setValueByGlobalID( gid[0], rhs );

        if ( fabs( pz - 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdz__( px, py, pz, 1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        } else if ( fabs( pz + 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdz__( px, py, pz, -1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        } else if ( fabs( px - 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdx__( px, py, pz, 1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        } else if ( fabs( px + 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdx__( px, py, pz, -1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        } else if ( fabs( py - 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdy__( px, py, pz, 1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        } else if ( fabs( py + 1.0 ) <= 1.0e-12 ) {
            val = __INIT_dTdy__( px, py, pz, -1.0 );
            val = val + __INIT_T0__( px, py, pz, 1.0 );
            variableFluxVec->setValueByGlobalID( gid[0], val );
        }
    } // end for node
    RightHandSideVec->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    correctionParameters->d_variableFlux = variableFluxVec;
    robinBoundaryOp->reset( correctionParameters );

    //----------------------------------------------------------
    //  Integrate Nuclear Rhs over Desnity * GeomType::Volume //
    //----------------------------------------------------------

    AMP_INSIST( input_db->keyExists( "VolumeIntegralOperator" ), "key missing!" );

    std::shared_ptr<AMP::Operator::ElementPhysicsModel> stransportModel;
    std::shared_ptr<AMP::Operator::VolumeIntegralOperator> sourceOperator =
        std::dynamic_pointer_cast<AMP::Operator::VolumeIntegralOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "VolumeIntegralOperator", input_db, stransportModel ) );

    // Create the power (heat source) vector.
    AMP::LinearAlgebra::Variable::shared_ptr SourceVar = sourceOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr SourceVec =
        AMP::LinearAlgebra::createVector( nodalDofMap, SourceVar, split );
    SourceVec->zero();

    // convert the vector of specific power to power for a given basis.
    sourceOperator->apply( RightHandSideVec, SourceVec );

    //------------------------------------------
    //   Add the boundary conditions corrections //
    //------------------------------------------

    std::cout << "RHS Norm before BC Correction " << SourceVec->L2Norm() << std::endl;

    diffusionOperator->modifyRHSvector( SourceVec );

    auto rhsNorm = SourceVec->L2Norm();
    std::cout << "RHS Norm after BC Correction " << rhsNorm << std::endl;

    //------------------------------------------

    // make sure the database on theinput file exists for the linear solver
    AMP_INSIST( input_db->keyExists( "LinearSolver" ), "Key ''LinearSolver'' is missing!" );

    // Read the input file onto a database.
    std::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "LinearSolver" );

    // Fill in the parameters fo the class with the info on the database.
    std::shared_ptr<AMP::Solver::SolverStrategyParameters> mlSolverParams(
        new AMP::Solver::SolverStrategyParameters( mlSolver_db ) );

    // Define the operature to be used by the Solver.
    mlSolverParams->d_pOperator = diffusionOperator;

    // Set initial guess
    TemperatureInKelvinVec->setToScalar( 1.0 );

    // Check the initial L2 norm of the solution
    double initSolNorm = TemperatureInKelvinVec->L2Norm();
    std::cout << "Initial Solution Norm: " << initSolNorm << std::endl;

    // Create the ML Solver
    std::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
        new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

    // Use a random initial guess?
    mlSolver->setZeroInitialGuess( false );

    // Solve the prblem.
    mlSolver->solve( SourceVec, TemperatureInKelvinVec );

    // Compute the residual
    diffusionOperator->residual( SourceVec, TemperatureInKelvinVec, ResidualVec );

    // Check the L2 norm of the final residual.
    double finalResidualNorm = ResidualVec->L2Norm();
    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    node                                            = node.begin();
    AMP::LinearAlgebra::Vector::shared_ptr diffVec  = TemperatureInKelvinVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr exactVec = TemperatureInKelvinVec->cloneVector();

    diffVec->zero();
    exactVec->zero();

    for ( ; node != end_node; ++node ) {
        std::vector<size_t> gid;
        nodalDofMap->getDOFs( node->globalID(), gid );

        double px = ( node->coord() )[0];
        double py = ( node->coord() )[1];
        double pz = ( node->coord() )[2];

        double exact;
        exact = __INIT_T0__( px, py, pz, 1.0 );
        exactVec->setValueByGlobalID( gid[0], exact );
    }

    diffVec->subtract( exactVec, TemperatureInKelvinVec );

    double exactNorm = exactVec->L1Norm();
    std::cout << "L2norm of exactVec " << exactNorm << std::endl;

    double solutionNorm = TemperatureInKelvinVec->L1Norm();
    std::cout << "L2norm of solutionVec " << solutionNorm << std::endl;

    double errorNorm = diffVec->L1Norm();
    std::cout << "L1norm of DiffVec " << errorNorm << std::endl;

    if ( errorNorm > 1.0 ) {
        ut->failure( "linear robin boundary operator verification test-1." );
    } else {
        ut->passes( "linear robin boundary operator verification test-1." );
    }

    // Plot the results
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );

#ifdef USE_EXT_SILO
    AMP::Utilities::Writer::shared_ptr siloWriter = AMP::Utilities::Writer::buildWriter( "Silo" );
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
