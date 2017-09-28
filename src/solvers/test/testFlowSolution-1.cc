#include "ampmesh/Mesh.h"
#include "operators/ElementOperationFactory.h"
#include "operators/ElementPhysicsModelFactory.h"
#include "operators/OperatorBuilder.h"
#include "operators/subchannel/FlowFrapconJacobian.h"
#include "operators/subchannel/FlowFrapconOperator.h"
#include "solvers/ColumnSolver.h"
#include "solvers/libmesh/Flow1DSolver.h"
#include "solvers/petsc/PetscKrylovSolver.h"
#include "solvers/petsc/PetscKrylovSolverParameters.h"
#include "solvers/petsc/PetscSNESSolver.h"
#include "solvers/petsc/PetscSNESSolverParameters.h"
#include "solvers/trilinos/ml/TrilinosMLSolver.h"
#include "utils/AMPManager.h"
#include "utils/Database.h"
#include "utils/InputDatabase.h"
#include "utils/InputManager.h"
#include "utils/PIO.h"
#include "utils/UnitTest.h"
#include "utils/Utilities.h"
#include "utils/shared_ptr.h"
#include "vectors/MultiVector.h"
#include "vectors/Variable.h"
#include "vectors/Vector.h"
#include <string>


void flowTest( AMP::UnitTest *ut, const std::string &exeName )
{
    std::string input_file = "input_" + exeName;
    std::string log_file   = "output_" + exeName;
    AMP::PIO::logAllNodes( log_file );
    AMP::AMP_MPI globalComm( AMP_COMM_WORLD );

    // Read the input file
    AMP::shared_ptr<AMP::InputDatabase> input_db( new AMP::InputDatabase( "input_db" ) );
    AMP::InputManager::getManager()->parseInputFile( input_file, input_db );

    // Get the Mesh database and create the mesh parameters
    AMP_INSIST( input_db->keyExists( "Mesh" ), "Key ''Mesh'' is missing!" );
    AMP::shared_ptr<AMP::Database> database = input_db->getDatabase( "Mesh" );
    AMP::shared_ptr<AMP::Mesh::MeshParameters> meshParams(
        new AMP::Mesh::MeshParameters( database ) );
    meshParams->setComm( globalComm );

    // Create the meshes from the input database
    AMP::shared_ptr<AMP::Mesh::Mesh> manager = AMP::Mesh::Mesh::buildMesh( meshParams );
    AMP::Mesh::Mesh::shared_ptr meshAdapter  = manager->Subset( "bar" );


    AMP::LinearAlgebra::Vector::shared_ptr nullVec;

    //--------------------------------------
    //     CREATE THE FLOW OPERATOR   ------
    //--------------------------------------

    AMP_INSIST( input_db->keyExists( "FlowFrapconOperator" ),
                "Key ''FlowFrapconOperator'' is missing!" );

    AMP::shared_ptr<AMP::Operator::ElementPhysicsModel> flowtransportModel;
    AMP::shared_ptr<AMP::InputDatabase> flowDatabase =
        AMP::dynamic_pointer_cast<AMP::InputDatabase>(
            input_db->getDatabase( "FlowFrapconOperator" ) );
    AMP::shared_ptr<AMP::Operator::FlowFrapconOperator> flowOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::FlowFrapconOperator>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "FlowFrapconOperator", input_db, flowtransportModel ) );

    AMP::LinearAlgebra::Variable::shared_ptr inputVariable  = flowOperator->getInputVariable();
    AMP::LinearAlgebra::Variable::shared_ptr outputVariable = flowOperator->getOutputVariable();

    AMP::LinearAlgebra::Vector::shared_ptr solVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr cladVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, inputVariable );

    AMP::LinearAlgebra::Vector::shared_ptr rhsVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr resVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, outputVariable );
    AMP::LinearAlgebra::Vector::shared_ptr tmpVec =
        AMP::LinearAlgebra::SimpleVector<double>::create( 10, inputVariable );

    AMP::shared_ptr<AMP::Operator::FlowFrapconJacobian> flowJacobian =
        AMP::dynamic_pointer_cast<AMP::Operator::FlowFrapconJacobian>(
            AMP::Operator::OperatorBuilder::createOperator(
                meshAdapter, "FlowFrapconJacobian", input_db, flowtransportModel ) );

    //-------------------------------------
    //    MANUFACTURE THE INPUT SOLUTION
    //     FROM PRESCRIBED FLOW SOLUTION

    // int cnt=0;
    double Tin = 300, Tc;
    double Cp, De, G, K, Re, Pr, heff, dz, nP;

    Cp = ( flowDatabase )->getDouble( "Heat_Capacity" );
    De = ( flowDatabase )->getDouble( "Channel_Diameter" );
    G  = ( flowDatabase )->getDouble( "Mass_Flux" );
    K  = ( flowDatabase )->getDouble( "Conductivity" );
    Re = ( flowDatabase )->getDouble( "Reynolds" );
    Pr = ( flowDatabase )->getDouble( "Prandtl" );
    nP = ( flowDatabase )->getDouble( "numpoints" );

    heff = ( 0.023 * K / De ) * pow( Re, 0.8 ) * pow( Pr, 0.4 );
    dz   = 2. / nP;

    std::cout << "Original Flow Solution  " << std::endl;
    for ( int i = 0; i < 10; i++ ) {
        tmpVec->setValueByLocalID( i, Tin + i * 30 );
        //	  solVec->setValueByLocalID(i, Tin + i*30);
        std::cout << " @i : " << i << " is " << tmpVec->getValueByLocalID( i );
    }
    std::cout << std::endl;

    std::cout << "Original Norm: " << tmpVec->L2Norm() << std::endl;

    cladVec->setValueByLocalID( 0, 300 );
    for ( int i = 1; i < 10; i++ ) {
        Tc = tmpVec->getValueByLocalID( i ) +
             ( tmpVec->getValueByLocalID( i ) - tmpVec->getValueByLocalID( i - 1 ) ) *
                 ( ( Cp * G * De ) / ( 4. * heff * dz ) );
        cladVec->setValueByLocalID( i, Tc );
    }

    std::cout << "Imposed Clad Solution  " << std::endl;
    for ( int i = 0; i < 10; i++ ) {
        std::cout << " @i : " << i << " is " << cladVec->getValueByLocalID( i );
    }
    std::cout << std::endl;

    flowOperator->setVector( cladVec );
    flowJacobian->setVector( cladVec );
    //  resVec->setToScalar(300.);
    //  vecLag->copyVector(resVec);
    //------------------------------------

    /*
      while ( 1 )
      {
          cnt++;
          flowOperator->apply(rhsVec, solVec, resVec, 1.0, -1.0);
    //	  if ( abs(resVec->L2Norm()-tmpVec->L2Norm()) < .000005*tmpVec->L2Norm() )
          if ( abs(resVec->L2Norm()-vecLag->L2Norm()) < .000005*vecLag->L2Norm() )
              break;
          else
              std::cout << "for iteration cnt = "<<cnt <<" --> " << vecLag->L2Norm() << " " <<
    resVec->L2Norm() <<
    std::endl;

          std::cout<<"Intermediate Flow Solution " <<std::endl;
          for( int i=0; i<10; i++) {
              std::cout<<" @i : "<< i<<" is "<<resVec->getValueByLocalID(i) ;
          }
          std::cout<<std::endl;
          vecLag->copyVector(resVec);
      }
    */
    AMP::shared_ptr<AMP::Database> jacobianSolver_db  = input_db->getDatabase( "JacobianSolver" );
    AMP::shared_ptr<AMP::Database> nonlinearSolver_db = input_db->getDatabase( "NonlinearSolver" );
    // AMP::shared_ptr<AMP::Database> linearSolver_db =
    // nonlinearSolver_db->getDatabase("LinearSolver");

    AMP::LinearAlgebra::Vector::shared_ptr mv_view_solVec =
        AMP::LinearAlgebra::MultiVector::view( solVec, globalComm );
    AMP::LinearAlgebra::Vector::shared_ptr mv_view_rhsVec =
        AMP::LinearAlgebra::MultiVector::view( rhsVec, globalComm );
    AMP::LinearAlgebra::Vector::shared_ptr mv_view_tmpVec =
        AMP::LinearAlgebra::MultiVector::view( tmpVec, globalComm );
    // AMP::LinearAlgebra::Vector::shared_ptr mv_view_cladVec =
    // AMP::LinearAlgebra::MultiVector::view( cladVec ,
    // globalComm );

    flowOperator->residual( rhsVec, solVec, resVec );
    flowJacobian->reset( flowOperator->getParameters( "Jacobian", mv_view_solVec ) );
    flowJacobian->residual( rhsVec, solVec, resVec );
    //----------------------------------------------------------------------------------------------------------------------------------------------//

    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> jacobianSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( jacobianSolver_db ) );

    // change the next line to get the correct communicator out
    jacobianSolverParams->d_comm          = globalComm;
    jacobianSolverParams->d_pOperator     = flowJacobian;
    jacobianSolverParams->d_pInitialGuess = mv_view_tmpVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> JacobianSolver(
        new AMP::Solver::PetscSNESSolver( jacobianSolverParams ) );

    //----------------------------------------------------------------------------------------------------------------------
    // initialize the nonlinear solver
    AMP::shared_ptr<AMP::Solver::PetscSNESSolverParameters> nonlinearSolverParams(
        new AMP::Solver::PetscSNESSolverParameters( nonlinearSolver_db ) );

    // change the next line to get the correct communicator out
    nonlinearSolverParams->d_comm          = globalComm;
    nonlinearSolverParams->d_pOperator     = flowOperator;
    nonlinearSolverParams->d_pInitialGuess = mv_view_tmpVec;

    AMP::shared_ptr<AMP::Solver::PetscSNESSolver> nonlinearSolver(
        new AMP::Solver::PetscSNESSolver( nonlinearSolverParams ) );

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    //  AMP::shared_ptr<AMP::Database> flowPreconditioner_db =
    //  linearSolver_db->getDatabase("Preconditioner");
    //  AMP::shared_ptr<AMP::Solver::SolverStrategyParameters> flowPreconditionerParams(new
    //  AMP::Solver::SolverStrategyParameters(flowPreconditioner_db));
    //  flowPreconditionerParams->d_pOperator = flowJacobian;
    //  AMP::shared_ptr<AMP::Solver::Flow1DSolver> linearFlowPreconditioner(new
    //  AMP::Solver::Flow1DSolver(flowPreconditionerParams));

    //----------------------------------------------------------------------------------------------------------------------------------------------//
    // register the preconditioner with the Jacobian free Krylov solver
    AMP::shared_ptr<AMP::Solver::PetscKrylovSolver> linearSolver =
        nonlinearSolver->getKrylovSolver();

    linearSolver->setPreconditioner( JacobianSolver );
    //------------------------------------------------------------------------------

    flowOperator->residual( rhsVec, solVec, resVec );
    double initialResidualNorm = resVec->L2Norm();

    AMP::pout << "Initial Residual Norm: " << initialResidualNorm << std::endl;

    nonlinearSolver->setZeroInitialGuess( true );

    nonlinearSolver->solve( mv_view_rhsVec, mv_view_solVec );

    flowOperator->residual( rhsVec, solVec, resVec );

    double finalResidualNorm = resVec->L2Norm();

    std::cout << "Final Residual Norm: " << finalResidualNorm << std::endl;

    std::cout << "Final Flow Solution " << std::endl;
    for ( int i = 0; i < 10; i++ ) {
        std::cout << " @i : " << i << " is " << solVec->getValueByLocalID( i );
    }
    std::cout << std::endl;
    //------------------------------------

    //  if( (tmpVec->L2Norm()-resVec->L2Norm()) > 0.01*tmpVec->L2Norm() )
    if ( resVec->L2Norm() > 0.01 ) {
        ut->failure( "Manufactured Solution verification test for 1D flow operator." );
    } else {
        ut->passes( "Manufactured Solution verification test for 1D flow operator." );
    }

    input_db.reset();

    ut->passes( exeName );
}


int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    flowTest( &ut, "testFlowSolution-1" );

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
