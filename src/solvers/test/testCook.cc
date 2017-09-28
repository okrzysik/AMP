#include "utils/AMPManager.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include <iostream>
#include <string>


#include <fstream>

#include <sys/stat.h>

#include "utils/shared_ptr.h"

// libMesh files
DISABLE_WARNINGS
#include "libmesh/boundary_info.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/elem.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_base.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_communication.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/quadrature.h"
#include "libmesh/string_to_enum.h"
ENABLE_WARNINGS

/* AMP files */
#include "utils/AMPManager.h"
#include "utils/AMP_MPI.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"

#include "ampmesh/Mesh.h"
#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "materials/Material.h"

#include "operators/ElementPhysicsModel.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"


#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "utils/Writer.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"
#include "vectors/VectorSelector.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

#include "utils/ReadTestMesh.h"

extern "C" {
#include "petsc.h"
}


void linearElasticTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName + ".txt";

    AMP::PIO::logOnlyNodeZero( log_file );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::AMP_MPI globalComm = AMP::AMP_MPI( AMP_COMM_WORLD );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );

    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( globalComm ) );

    const unsigned int mesh_dim = 3;
    AMP::shared_ptr<::Mesh> mesh( new ::Mesh( mesh_dim ) );

    std::string mesh_file = input_db->getString( "mesh_file" );
    if ( globalComm.getRank() == 0 ) {
        AMP::readTestMesh( mesh_file, mesh );
    } // end if root processor

    MeshCommunication().broadcast( *( mesh.get() ) );
    mesh->prepare_for_use( false );
    AMP::Mesh::Mesh::shared_ptr meshAdapter( new AMP::Mesh::libMesh( mesh, "cook" ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> elementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> bvpOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "MechanicsBVPOperator", input_db, elementPhysicsModel ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> dirichletVecOp =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "Load_Boundary", input_db, dummyModel ) );
    // This has an in-place apply. So, it has an empty input variable and
    // the output variable is the same as what it is operating on.
    dirichletVecOp->setVariable( bvpOperator->getOutputVariable() );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    AMP::Discretization::DOFManager::shared_ptr DOF_vector =
        AMP::Discretization::simpleDOFManager::create(
            meshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3, true );
    AMP::LinearAlgebra::Vector::shared_ptr mechSolVec =
        AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );
    AMP::LinearAlgebra::Vector::shared_ptr mechRhsVec =
        AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );
    AMP::LinearAlgebra::Vector::shared_ptr mechResVec =
        AMP::LinearAlgebra::createVector( DOF_vector, bvpOperator->getOutputVariable(), true );

    mechSolVec->setToScalar( 0.5 );
    mechRhsVec->setToScalar( 0.0 );
    mechResVec->setToScalar( 0.0 );

    dirichletVecOp->apply( nullVec, mechRhsVec );

    double rhsNorm = mechRhsVec->L2Norm();

    AMP::pout << "RHS Norm: " << rhsNorm << std::endl;

    double initSolNorm = mechSolVec->L2Norm();

    AMP::pout << "Initial Solution Norm: " << initSolNorm << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double initResidualNorm = mechResVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initResidualNorm << std::endl;

    AMP::shared_ptr<AMP::Database> linearSolver_db = input_db->getDatabase( "LinearSolver" );

    // ---- first initialize the preconditioner
    AMP::shared_ptr<AMP::Database> pcSolver_db = linearSolver_db->getDatabase( "Preconditioner" );
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> pcSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( pcSolver_db ) );
    pcSolverParams->d_pOperator = bvpOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> pcSolver(
        new AMP::Solver::TrilinosMLSolver( pcSolverParams ) );

    // initialize the linear solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolverParameters> linearSolverParams(
        new AMP::Solver::PetscKrylovSolverParameters( linearSolver_db ) );
    linearSolverParams->d_pOperator       = bvpOperator;
    linearSolverParams->d_comm            = globalComm;
    linearSolverParams->d_pPreconditioner = pcSolver;
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver(
        new AMP::Solver::PetscKrylovSolver( linearSolverParams ) );

    linearSolver->setZeroInitialGuess( false );

    linearSolver->solve( mechRhsVec, mechSolVec );

    double finalSolNorm = mechSolVec->L2Norm();

    AMP::pout << "Final Solution Norm: " << finalSolNorm << std::endl;

    AMP::LinearAlgebra::Vector::shared_ptr mechUvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 0, 3 ), "U" );
    AMP::LinearAlgebra::Vector::shared_ptr mechVvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 1, 3 ), "V" );
    AMP::LinearAlgebra::Vector::shared_ptr mechWvec =
        mechSolVec->select( AMP::LinearAlgebra::VS_Stride( 2, 3 ), "W" );

    double finalMaxU = mechUvec->maxNorm();
    double finalMaxV = mechVvec->maxNorm();
    double finalMaxW = mechWvec->maxNorm();

    AMP::pout << "Maximum U displacement: " << finalMaxU << std::endl;
    AMP::pout << "Maximum V displacement: " << finalMaxV << std::endl;
    AMP::pout << "Maximum W displacement: " << finalMaxW << std::endl;

    bvpOperator->residual( mechRhsVec, mechSolVec, mechResVec );

    double finalResidualNorm = mechResVec->L2Norm();

    AMP::pout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    if ( finalResidualNorm > ( 1e-10 * initResidualNorm ) ) {
        ut->failure( exeName );
    } else {
        ut->passes( exeName );
    }
}

int main( int argc, char *argv[] )
{

    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "testCook-normal-mesh0" );
        exeNames.emplace_back( "testCook-reduced-mesh0" );

        exeNames.emplace_back( "testCook-normal-mesh1" );
        exeNames.emplace_back( "testCook-reduced-mesh1" );

        exeNames.emplace_back( "testCook-normal-mesh2" );
        exeNames.emplace_back( "testCook-reduced-mesh2" );
    } else {
        for ( int i = 1; i < argc; i += 2 ) {
            char inpName[100];
            sprintf( inpName, "testCook-%s-mesh%d", argv[i], atoi( argv[i + 1] ) );
            exeNames.emplace_back( inpName );
        } // end for i
    }

    for ( auto &exeName : exeNames )
        linearElasticTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
