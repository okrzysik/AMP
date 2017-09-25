#include <iostream>
#include <string>

#include "materials/Material.h"

#include "utils/AMPManager.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/ReadTestMesh.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/WriteSolutionToFile.h"

#include "ampmesh/MultiMesh.h"
#include "ampmesh/libmesh/initializeLibMesh.h"
#include "ampmesh/libmesh/libMesh.h"
#include "utils/ReadTestMesh.h"

#include "discretization/DOF_Manager.h"
#include "discretization/simpleDOF_Manager.h"
#include "vectors/Vector.h"
#include "vectors/VectorBuilder.h"

#include "libmesh/mesh_communication.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearBVPOperator.h"
#include "operators/OperatorBuilder.h"
#include "operators/boundary/DirichletVectorCorrection.h"
#include "operators/trilinos/TrilinosMatrixShellOperator.h"

#include "vectors/trilinos/epetra/EpetraVector.h"

#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"

#include "ml_include.h"

void myGetRow2( void *object, int row, std::vector<size_t> &cols, std::vector<double> &values )
{
    auto *op = reinterpret_cast<AMP::Operator::ColumnOperator *>( object );
    AMP::LinearAlgebra::Matrix::shared_ptr mat =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op->getOperator( 0 ) )
            ->getMatrix();
    mat->getRowByGlobalID( row, cols, values );
}

void myGetRow3( void *object, int row, std::vector<size_t> &cols, std::vector<double> &values )
{
    auto *op = reinterpret_cast<AMP::Operator::ColumnOperator *>( object );
    AMP::LinearAlgebra::Matrix::shared_ptr firstMat =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op->getOperator( 0 ) )
            ->getMatrix();
    AMP::LinearAlgebra::Matrix::shared_ptr secondMat =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op->getOperator( 1 ) )
            ->getMatrix();
    size_t firstMatNumGlobalRows    = firstMat->numGlobalRows();
    size_t firstMatNumGlobalColumns = firstMat->numGlobalColumns();
    if ( row < (int) firstMatNumGlobalRows ) {
        firstMat->getRowByGlobalID( row, cols, values );
    } else {
        secondMat->getRowByGlobalID( row - firstMatNumGlobalRows, cols, values );
        for ( auto &col : cols ) {
            col += firstMatNumGlobalColumns;
        } // end for j
    }     // end if
}

int myMatVec( ML_Operator *data, int in_length, double in[], int out_length, double out[] )
{

    auto *op = reinterpret_cast<AMP::Operator::LinearOperator *>( ML_Get_MyMatvecData( data ) );
    AMP::LinearAlgebra::Matrix::shared_ptr mat = op->getMatrix();

    AMP::LinearAlgebra::Vector::shared_ptr inVec  = mat->getRightVector();
    AMP::LinearAlgebra::Vector::shared_ptr outVec = mat->getLeftVector();

    AMP_ASSERT( in_length == (int) inVec->getLocalSize() );
    AMP_ASSERT( out_length == (int) outVec->getLocalSize() );

    inVec->putRawData( in );

    mat->mult( inVec, outVec );

    outVec->copyOutRawData( out );

    return 0;
}


int myGetRow( ML_Operator *data,
              int N_requested_rows,
              int requested_rows[],
              int allocated_space,
              int columns[],
              double values[],
              int row_lengths[] )
{

    auto *op = reinterpret_cast<AMP::Operator::LinearOperator *>( ML_Get_MyGetrowData( data ) );
    AMP::LinearAlgebra::Matrix::shared_ptr mat = op->getMatrix();

    int spaceRequired = 0;
    int cnt           = 0;
    for ( int i = 0; i < N_requested_rows; i++ ) {
        int row = requested_rows[i];
        std::vector<size_t> cols;
        std::vector<double> vals;

        mat->getRowByGlobalID( row, cols, vals );
        spaceRequired += cols.size();

        if ( allocated_space >= spaceRequired ) {
            for ( size_t j = 0; j < cols.size(); j++ ) {
                columns[cnt] = cols[j];
                values[cnt]  = vals[j];
                cnt++;
            }
            row_lengths[i] = cols.size();
        } else {
            return 0;
        }
    }

    return 1;
}


void myTest( AMP::UnitTest *ut, std::string exeName, int type )
{
    std::string input_file = "input_" + exeName;
    char log_file[200];
    sprintf( log_file, "output_%s_%d", exeName.c_str(), type );

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );
    std::string mesh_file = input_db->getString( "mesh_file" );

    AMP::shared_ptr<AMP::InputDatabase> mesh_file_db( new AMP::InputDatabase( "mesh_file_db" ) );
    AMP::InputManager::getManager()->parseInputFile( mesh_file, mesh_file_db );

    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( globalComm ) );

    const unsigned int mesh_dim = 3;
    AMP::shared_ptr<::Mesh> fusedMesh( new ::Mesh( mesh_dim ) );

    AMP::readTestMesh( mesh_file, fusedMesh );

    MeshCommunication().broadcast( *( fusedMesh.get() ) );

    fusedMesh->prepare_for_use( false );

    AMP::Mesh::Mesh::shared_ptr fusedMeshAdapter( new AMP::Mesh::libMesh( fusedMesh, "mesh" ) );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> fusedElementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> fusedOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                fusedMeshAdapter, "BVPOperator", input_db, fusedElementPhysicsModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr fusedVar = fusedOperator->getOutputVariable();

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> loadOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                fusedMeshAdapter, "LoadOperator", input_db, dummyModel ) );
    loadOperator->setVariable( fusedVar );

    AMP::Discretization::DOFManager::shared_ptr NodalVectorDOF =
        AMP::Discretization::simpleDOFManager::create(
            fusedMeshAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );

    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr fusedSolVec =
        AMP::LinearAlgebra::createVector( NodalVectorDOF, fusedVar );
    AMP::LinearAlgebra::Vector::shared_ptr fusedRhsVec = fusedSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr fusedResVec = fusedSolVec->cloneVector();

    fusedRhsVec->zero();
    loadOperator->apply( nullVec, fusedRhsVec );

    AMP::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "MLoptions" );

    std::cout << std::endl;

    const size_t localSize = fusedSolVec->getLocalSize();
    NULL_USE( localSize );

    // Matrix-based
    if ( type == 0 ) {
        ML_set_random_seed( 123456 );
        std::cout << "Matrix-Based ML: " << std::endl;
        fusedSolVec->zero();

        AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlSolverParams(
            new AMP::Solver::TrilinosMLSolverParameters( mlSolver_db ) );
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
            new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

        AMP::LinearAlgebra::Matrix::shared_ptr mat = fusedOperator->getMatrix();

        AMP::LinearAlgebra::Matrix::shared_ptr matCopy = mat->cloneMatrix();
        matCopy->zero();
        matCopy->axpy( 1.0, mat );

        mat->zero();
        mlSolver->registerOperator( fusedOperator );

        mat->axpy( 1.0, matCopy );

        mlSolver->solve( fusedRhsVec, fusedSolVec );
        std::cout << std::endl;
    }

    // Matrix-Free-1
    if ( type == 1 ) {
        ML_set_random_seed( 123456 );
        std::cout << "Matrix-Free ML Type-1: " << std::endl;
        fusedSolVec->zero();

        ML_Comm *comm;
        ML_Comm_Create( &comm );
        ML_Comm_Set_UsrComm( comm, globalComm.getCommunicator() );

        ML_Operator *ml_op = ML_Operator_Create( comm );
        ML_Operator_Set_ApplyFuncData(
            ml_op, localSize, localSize, fusedOperator.get(), localSize, myMatVec, 0 );
        ML_Operator_Set_Getrow( ml_op, localSize, myGetRow );

        Teuchos::ParameterList paramsList;
        ML_Epetra::SetDefaults( "SA", paramsList );
        paramsList.set( "ML output", mlSolver_db->getInteger( "print_info_level" ) );
        paramsList.set( "PDE equations", mlSolver_db->getInteger( "PDE_equations" ) );
        paramsList.set( "cycle applications", mlSolver_db->getInteger( "max_iterations" ) );
        paramsList.set( "max levels", mlSolver_db->getInteger( "max_levels" ) );

        AMP::shared_ptr<ML_Epetra::MultiLevelPreconditioner> mlSolver(
            new ML_Epetra::MultiLevelPreconditioner( ml_op, paramsList ) );

        const ML_Aggregate *agg_obj = mlSolver->GetML_Aggregate();
        ML_Aggregate_Print( const_cast<ML_Aggregate *>( agg_obj ) );

        auto f_epetra = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraVector>(
            AMP::LinearAlgebra::EpetraVector::view( fusedRhsVec ) );
        auto u_epetra = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::EpetraVector>(
            AMP::LinearAlgebra::EpetraVector::view( fusedSolVec ) );

        Epetra_Vector &fVec = f_epetra->getEpetra_Vector();
        Epetra_Vector &uVec = u_epetra->getEpetra_Vector();

        fusedOperator->residual( fusedRhsVec, fusedSolVec, fusedResVec );
        std::cout << "MatFree-1: L2 norm of residual before solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        mlSolver->ApplyInverse( fVec, uVec );

        auto firer = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::DataChangeFirer>( fusedSolVec );
        if ( firer )
            firer->fireDataChange();

        double solution_norm = fusedSolVec->L2Norm();
        std::cout << "MatFree-1:  solution norm: " << std::setprecision( 15 ) << solution_norm
                  << std::endl;

        fusedOperator->residual( fusedRhsVec, fusedSolVec, fusedResVec );
        std::cout << "MatFree-1: L2 norm of residual after solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        ML_Operator_Destroy( &ml_op );

        ML_Comm_Destroy( &comm );
    }

    // Matrix-free-2
    if ( type == 2 ) {
        ML_set_random_seed( 123456 );
        std::cout << "Matrix-Free ML Type-2: " << std::endl;
        fusedSolVec->zero();

        int numGrids = mlSolver_db->getInteger( "max_levels" );
        int numPDEs  = mlSolver_db->getInteger( "PDE_equations" );

        ML *ml_object;
        ML_Create( &ml_object, numGrids );

        ML_Init_Amatrix( ml_object, 0, localSize, localSize, fusedOperator.get() );
        ML_Set_Amatrix_Getrow( ml_object, 0, &myGetRow, nullptr, localSize );
        ML_Set_Amatrix_Matvec( ml_object, 0, &myMatVec );
        ML_Set_MaxIterations( ml_object, 1 + mlSolver_db->getInteger( "max_iterations" ) );
        ML_Set_ResidualOutputFrequency( ml_object, 1 );
        ML_Set_PrintLevel( mlSolver_db->getInteger( "print_info_level" ) );
        ML_Set_OutputLevel( ml_object, mlSolver_db->getInteger( "print_info_level" ) );

        ML_Aggregate *agg_object;
        ML_Aggregate_Create( &agg_object );
        agg_object->num_PDE_eqns  = numPDEs;
        agg_object->nullspace_dim = numPDEs;
        ML_Aggregate_Set_MaxCoarseSize( agg_object, 128 );
        ML_Aggregate_Set_CoarsenScheme_UncoupledMIS( agg_object );

        int nlevels =
            ML_Gen_MGHierarchy_UsingAggregation( ml_object, 0, ML_INCREASING, agg_object );
        std::cout << "Number of actual levels : " << nlevels << std::endl;

        for ( int lev = 0; lev < ( nlevels - 1 ); lev++ ) {
            ML_Gen_Smoother_SymGaussSeidel( ml_object, lev, ML_BOTH, 2, 1.0 );
        }
        ML_Gen_Smoother_Amesos( ml_object, ( nlevels - 1 ), ML_AMESOS_KLU, -1, 0.0 );

        ML_Gen_Solver( ml_object, ML_MGV, 0, ( nlevels - 1 ) );


        fusedOperator->residual( fusedRhsVec, fusedSolVec, fusedResVec );
        std::cout << "MatFree-2: L2 norm of residual before solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        auto solArr = new double[fusedSolVec->getLocalSize()];
        auto rhsArr = new double[fusedRhsVec->getLocalSize()];
        fusedSolVec->copyOutRawData( solArr );
        fusedRhsVec->copyOutRawData( rhsArr );
        ML_Iterate( ml_object, solArr, rhsArr );
        fusedSolVec->putRawData( solArr );
        fusedRhsVec->putRawData( rhsArr );
        delete[] solArr;
        solArr = nullptr;
        delete[] rhsArr;
        rhsArr = nullptr;

        auto firer = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::DataChangeFirer>( fusedSolVec );
        if ( firer )
            firer->fireDataChange();

        double solution_norm = fusedSolVec->L2Norm();
        std::cout << "MatFree-2:  solution norm: " << std::setprecision( 15 ) << solution_norm
                  << std::endl;

        fusedOperator->residual( fusedRhsVec, fusedSolVec, fusedResVec );
        std::cout << "MatFree-2: L2 norm of residual after solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        ML_Aggregate_Destroy( &agg_object );
        ML_Destroy( &ml_object );
        std::cout << std::endl;
    }

    // matrix-free-3 using TrilinosMatrixShellOperator and customized getRow()
    if ( type == 3 ) {
        AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyParams;
        AMP::shared_ptr<AMP::Operator::ColumnOperator> columnOperator(
            new AMP::Operator::ColumnOperator( emptyParams ) );
        columnOperator->append( fusedOperator );

        AMP::shared_ptr<AMP::Database> matrixShellDatabase(
            new AMP::MemoryDatabase( "MatrixShellOperator" ) );
        matrixShellDatabase->putString( "name", "MatShellOperator" );
        matrixShellDatabase->putInteger( "print_info_level", 1 );
        AMP::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(
            new AMP::Operator::OperatorParameters( matrixShellDatabase ) );
        AMP::shared_ptr<AMP::Operator::TrilinosMatrixShellOperator> trilinosMatrixShellOperator(
            new AMP::Operator::TrilinosMatrixShellOperator( matrixShellParams ) );
        trilinosMatrixShellOperator->setNodalDofMap( NodalVectorDOF );
        trilinosMatrixShellOperator->setGetRow( &myGetRow2 );
        // trilinosMatrixShellOperator->setOperator(fusedOperator);
        trilinosMatrixShellOperator->setOperator( columnOperator );

        mlSolver_db->putBool( "USE_EPETRA", false );
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlSolverParams(
            new AMP::Solver::TrilinosMLSolverParameters( mlSolver_db ) );
        mlSolverParams->d_pOperator = trilinosMatrixShellOperator;
        AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
            new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

        std::cout << "MatFree-3: L2 norm of residual before solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        mlSolver->solve( fusedRhsVec, fusedSolVec );

        std::cout << "MatFree-3:  solution norm: " << std::setprecision( 15 )
                  << fusedSolVec->L2Norm() << std::endl;
        fusedOperator->residual( fusedRhsVec, fusedSolVec, fusedResVec );
        std::cout << "MatFree-3: L2 norm of residual after solve " << std::setprecision( 15 )
                  << fusedResVec->L2Norm() << std::endl;

        std::cout << std::endl;
    }

    char outFile[200];
    sprintf( outFile, "%s-%d", exeName.c_str(), type );
    printSolution( fusedMeshAdapter, fusedSolVec, outFile );

    ut->passes( exeName );
}

void myTest2( AMP::UnitTest *ut, std::string exeName, bool useTwoMeshes )
{
    std::string input_file = "input_" + exeName;
    char log_file[200];
    int type = 4;
    if ( useTwoMeshes ) {
        type = 5;
    } // end if
    sprintf( log_file, "output_%s_%d", exeName.c_str(), type );

    AMP::PIO::logOnlyNodeZero( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );
    input_db->printClassData( AMP::plog );
    std::string mesh_file = input_db->getString( "mesh_file" );

    AMP::shared_ptr<AMP::InputDatabase> mesh_file_db( new AMP::InputDatabase( "mesh_file_db" ) );
    AMP::InputManager::getManager()->parseInputFile( mesh_file, mesh_file_db );

    AMP::shared_ptr<AMP::Mesh::initializeLibMesh> libmeshInit(
        new AMP::Mesh::initializeLibMesh( globalComm ) );

    const unsigned int mesh_dim = 3;
    AMP::shared_ptr<::Mesh> firstFusedMesh( new ::Mesh( mesh_dim ) );
    AMP::shared_ptr<::Mesh> secondFusedMesh( new ::Mesh( mesh_dim ) );

    AMP::readTestMesh( mesh_file, firstFusedMesh );
    AMP::readTestMesh( mesh_file, secondFusedMesh );

    MeshCommunication().broadcast( *( firstFusedMesh.get() ) );
    MeshCommunication().broadcast( *( secondFusedMesh.get() ) );

    firstFusedMesh->prepare_for_use( false );
    secondFusedMesh->prepare_for_use( false );

    AMP::Mesh::Mesh::shared_ptr firstMesh( new AMP::Mesh::libMesh( firstFusedMesh, "Mesh_1" ) );
    AMP::Mesh::Mesh::shared_ptr secondMesh( new AMP::Mesh::libMesh( secondFusedMesh, "Mesh_2" ) );

    std::vector<AMP::Mesh::Mesh::shared_ptr> vectorOfMeshes;
    vectorOfMeshes.push_back( firstMesh );
    if ( useTwoMeshes ) {
        vectorOfMeshes.push_back( secondMesh );
    } // end if

    AMP::Mesh::Mesh::shared_ptr fusedMeshesAdapter(
        new AMP::Mesh::MultiMesh( globalComm, vectorOfMeshes ) );
    std::vector<AMP::Mesh::Mesh::shared_ptr> fusedMeshes =
        AMP::dynamic_pointer_cast<AMP::Mesh::MultiMesh>( fusedMeshesAdapter )->getMeshes();
    fusedMeshesAdapter->setName( "MultiMesh" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> fusedElementPhysicsModel;
    AMP::shared_ptr<AMP::Operator::LinearBVPOperator> firstFusedOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                fusedMeshes[0], "BVPOperator", input_db, fusedElementPhysicsModel ) );
    AMP::LinearAlgebra::Variable::shared_ptr firstFusedVar =
        firstFusedOperator->getOutputVariable();
    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> dummyModel;
    AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> firstLoadOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
            AMP::Operator::OperatorBuilder::createOperator(
                fusedMeshes[0], "LoadOperator", input_db, dummyModel ) );
    firstLoadOperator->setVariable( firstFusedVar );

    AMP::shared_ptr<AMP::Operator::OperatorParameters> emptyOperatorParams;
    AMP::shared_ptr<AMP::Operator::ColumnOperator> fusedColumnOperator(
        new AMP::Operator::ColumnOperator( emptyOperatorParams ) );
    fusedColumnOperator->append( firstFusedOperator );

    AMP::shared_ptr<AMP::Operator::ColumnOperator> loadColumnOperator(
        new AMP::Operator::ColumnOperator( emptyOperatorParams ) );
    loadColumnOperator->append( firstLoadOperator );
    if ( useTwoMeshes ) {
        AMP::shared_ptr<AMP::Operator::LinearBVPOperator> secondFusedOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::LinearBVPOperator>(
                AMP::Operator::OperatorBuilder::createOperator(
                    fusedMeshes[1], "BVPOperator", input_db, fusedElementPhysicsModel ) );
        AMP::LinearAlgebra::Variable::shared_ptr secondFusedVar =
            secondFusedOperator->getOutputVariable();
        AMP::shared_ptr<AMP::Operator::DirichletVectorCorrection> secondLoadOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::DirichletVectorCorrection>(
                AMP::Operator::OperatorBuilder::createOperator(
                    fusedMeshes[1], "LoadOperator", input_db, dummyModel ) );
        secondLoadOperator->setVariable( secondFusedVar );
        fusedColumnOperator->append( secondFusedOperator );
        loadColumnOperator->append( secondLoadOperator );
    } // end if


    AMP::Discretization::DOFManager::shared_ptr dofManager =
        AMP::Discretization::simpleDOFManager::create(
            fusedMeshesAdapter, AMP::Mesh::GeomType::Vertex, 1, 3 );

    AMP::LinearAlgebra::Variable::shared_ptr fusedColumnVar =
        fusedColumnOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr nullVec;
    AMP::LinearAlgebra::Vector::shared_ptr fusedColumnSolVec =
        AMP::LinearAlgebra::createVector( dofManager, fusedColumnVar );
    AMP::LinearAlgebra::Vector::shared_ptr fusedColumnRhsVec = fusedColumnSolVec->cloneVector();
    AMP::LinearAlgebra::Vector::shared_ptr fusedColumnResVec = fusedColumnSolVec->cloneVector();

    fusedColumnRhsVec->zero();
    loadColumnOperator->apply( nullVec, fusedColumnRhsVec );

    std::cout << std::endl;

    AMP::shared_ptr<AMP::Database> matrixShellDatabase(
        new AMP::MemoryDatabase( "MatrixShellOperator" ) );
    matrixShellDatabase->putString( "name", "MatShellOperator" );
    matrixShellDatabase->putInteger( "print_info_level", 1 );
    AMP::shared_ptr<AMP::Operator::OperatorParameters> matrixShellParams(
        new AMP::Operator::OperatorParameters( matrixShellDatabase ) );
    AMP::shared_ptr<AMP::Operator::TrilinosMatrixShellOperator> trilinosMatrixShellOperator(
        new AMP::Operator::TrilinosMatrixShellOperator( matrixShellParams ) );
    trilinosMatrixShellOperator->setNodalDofMap( dofManager );
    if ( !useTwoMeshes ) {
        trilinosMatrixShellOperator->setGetRow( &myGetRow2 );
    } else {
        trilinosMatrixShellOperator->setGetRow( &myGetRow3 );
    } // end if
    trilinosMatrixShellOperator->setOperator( fusedColumnOperator );


    AMP::shared_ptr<AMP::Database> mlSolver_db = input_db->getDatabase( "MLoptions" );
    mlSolver_db->putBool( "USE_EPETRA", false );
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolverParameters> mlSolverParams(
        new AMP::Solver::TrilinosMLSolverParameters( mlSolver_db ) );
    mlSolverParams->d_pOperator = trilinosMatrixShellOperator;
    AMP::shared_ptr<AMP::Solver::TrilinosMLSolver> mlSolver(
        new AMP::Solver::TrilinosMLSolver( mlSolverParams ) );

    std::cout << "MatFree-4: L2 norm of residual before solve " << std::setprecision( 15 )
              << fusedColumnResVec->L2Norm() << std::endl;

    mlSolver->solve( fusedColumnRhsVec, fusedColumnSolVec );

    std::cout << "MatFree-4:  solution norm: " << std::setprecision( 15 )
              << fusedColumnSolVec->L2Norm() << std::endl;
    fusedColumnOperator->residual( fusedColumnRhsVec, fusedColumnSolVec, fusedColumnResVec );
    std::cout << "MatFree-4: L2 norm of residual after solve " << std::setprecision( 15 )
              << fusedColumnResVec->L2Norm() << std::endl;

    std::cout << std::endl;

    ut->passes( exeName );
}


void loopMyTest( AMP::UnitTest *ut, std::string exeName )
{
    myTest2( ut, exeName, false );
    myTest2( ut, exeName, true );
    for ( int type = 0; type < 4; type++ ) {
        myTest( ut, exeName, type );
    }
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    std::vector<std::string> exeNames;

    if ( argc == 1 ) {
        exeNames.emplace_back( "testMatrixFreeML-1" );
    } else {
        for ( int i = 1; i < argc; i++ ) {
            char inpName[100];
            sprintf( inpName, "testMatrixFreeML-%s", argv[i] );
            exeNames.emplace_back( inpName );
        } // end for i
    }

    for ( auto &exeName : exeNames )
        loopMyTest( &ut, exeName );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
