#include "AMP/solvers/trilinos/nox/TrilinosNOXSolver.h"
#include "ProfilerApp.h"

#include "AMP/solvers/trilinos/nox/AndersonStatusTest.h"
#include "AMP/solvers/trilinos/thyra/TrilinosThyraModelEvaluator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/trilinos/thyra/ThyraVector.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


// Trilinos includes
DISABLE_WARNINGS

#ifdef USE_TRILINOS_BELOS
#include "BelosTypes.hpp"
#endif

#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"
#include "NOX_Solver_Factory.H"
#include "NOX_StatusTest_Combo.H"
#include "NOX_StatusTest_FiniteValue.H"
#include "NOX_StatusTest_MaxIters.H"
#include "NOX_StatusTest_NormF.H"
#include "NOX_StatusTest_NormUpdate.H"
#include "NOX_StatusTest_NormWRMS.H"
#include "NOX_StatusTest_RelativeNormF.H"
#include "NOX_Thyra.H"
#include "NOX_Thyra_Group.H"
#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include <Teuchos_RefCountPtrDecl.hpp>
ENABLE_WARNINGS


namespace AMP {
namespace Solver {


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
TrilinosNOXSolver::TrilinosNOXSolver() : SolverStrategy() {}
TrilinosNOXSolver::TrilinosNOXSolver( AMP::shared_ptr<TrilinosNOXSolverParameters> parameters )
    : SolverStrategy( parameters )
{
    initialize( parameters );
}
void TrilinosNOXSolver::reset( AMP::shared_ptr<SolverStrategyParameters> parameters )
{
    initialize( parameters );
}
TrilinosNOXSolver::~TrilinosNOXSolver() = default;


/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void TrilinosNOXSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> parameters )
{
    // Copy the parameters
    AMP::shared_ptr<TrilinosNOXSolverParameters> params =
        AMP::dynamic_pointer_cast<TrilinosNOXSolverParameters>( parameters );
    AMP_ASSERT( params.get() != nullptr );
    AMP_ASSERT( params->d_db.get() != nullptr );
    d_comm = params->d_comm;
    if ( params->d_pInitialGuess.get() != nullptr )
        d_initialGuess = params->d_pInitialGuess;
    AMP_ASSERT( d_initialGuess != nullptr );
    AMP::shared_ptr<AMP::Database> nonlinear_db = parameters->d_db;
    AMP::shared_ptr<AMP::Database> linear_db    = nonlinear_db->getDatabase( "LinearSolver" );
    AMP_ASSERT( linear_db != nullptr );
    // Create a model evaluator
    auto modelParams           = AMP::make_shared<TrilinosThyraModelEvaluatorParameters>();
    modelParams->d_nonlinearOp = d_pOperator;
    modelParams->d_linearOp    = params->d_pLinearOperator;
    modelParams->d_icVec       = d_initialGuess;
    modelParams->d_preconditioner.reset();
    modelParams->d_prePostOperator = params->d_prePostOperator;
    if ( linear_db->getBoolWithDefault( "uses_preconditioner", false ) )
        modelParams->d_preconditioner = params->d_preconditioner;
    d_thyraModel =
        Teuchos::RCP<TrilinosThyraModelEvaluator>( new TrilinosThyraModelEvaluator( modelParams ) );
    // Create the Preconditioner operator
    d_precOp = d_thyraModel->create_W_prec();
    // Create the linear solver factory
    ::Stratimikos::DefaultLinearSolverBuilder builder;
    Teuchos::RCP<Teuchos::ParameterList> p( new Teuchos::ParameterList );
    std::string linearSolverType   = linear_db->getString( "linearSolverType" );
    std::string linearSolver       = linear_db->getString( "linearSolver" );
    int maxLinearIterations        = linear_db->getIntegerWithDefault( "max_iterations", 100 );
    double linearRelativeTolerance = linear_db->getDoubleWithDefault( "relative_tolerance", 1e-3 );
    bool flexGmres                 = linear_db->getBoolWithDefault( "flexibleGmres", true );
    p->set( "Linear Solver Type", linearSolverType );
    p->set( "Preconditioner Type", "None" );
    p->sublist( "Linear Solver Types" )
        .sublist( linearSolverType )
        .set( "Solver Type", linearSolver );
    Teuchos::ParameterList &linearSolverParams =
        p->sublist( "Linear Solver Types" ).sublist( linearSolverType );
    linearSolverParams.sublist( "Solver Types" )
        .sublist( linearSolver )
        .set( "Maximum Iterations", maxLinearIterations );
    // Only "Block GMRES" recognizes the "Flexible Gmres" option, other solvers may throw an input
    // validation error
    if ( linearSolver == "Block GMRES" )
        linearSolverParams.sublist( "Solver Types" )
            .sublist( linearSolver )
            .set( "Flexible Gmres", flexGmres );
    if ( linear_db->getIntegerWithDefault( "print_info_level", 0 ) >= 2 ) {
        linearSolverParams.sublist( "Solver Types" )
            .sublist( linearSolver )
            .set( "Output Frequency", 1 );
        linearSolverParams.sublist( "Solver Types" ).sublist( linearSolver ).set( "Verbosity", 10 );
        linearSolverParams.sublist( "VerboseObject" ).set( "Verbosity Level", "extreme" );
        if ( linearSolverType == "Belos" ) {
            linearSolverParams.sublist( "Solver Types" )
                .sublist( linearSolver )
                .set( "Verbosity",
                      Belos::Warnings + Belos::IterationDetails + Belos::OrthoDetails +
                          Belos::FinalSummary + Belos::Debug + Belos::StatusTestDetails );
        }
    } else if ( linear_db->getIntegerWithDefault( "print_info_level", 0 ) >= 1 ) {
        linearSolverParams.sublist( "Solver Types" )
            .sublist( linearSolver )
            .set( "Output Frequency", 1 );
        linearSolverParams.sublist( "Solver Types" ).sublist( linearSolver ).set( "Verbosity", 10 );
        if ( linearSolverType == "Belos" ) {
            linearSolverParams.sublist( "Solver Types" )
                .sublist( linearSolver )
                .set( "Verbosity",
                      Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary +
                          Belos::Debug );
        }
    }
    builder.setParameterList( p );
    d_lowsFactory = builder.createLinearSolveStrategy( "" );
    // d_lowsFactory->initializeVerboseObjectBase();
    d_thyraModel->set_W_factory( d_lowsFactory );
    // Create the convergence tests (these will need to be on the input database)
    Teuchos::RCP<NOX::StatusTest::NormF> absresid( new NOX::StatusTest::NormF( d_dMaxError ) );
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters(
        new NOX::StatusTest::MaxIters( d_iMaxIterations ) );
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv( new NOX::StatusTest::FiniteValue );
    Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms(
        new NOX::StatusTest::NormWRMS( d_dMaxError, d_dMaxError ) );
    d_status = Teuchos::rcp( new NOX::StatusTest::Combo( NOX::StatusTest::Combo::OR ) );
    d_status->addStatusTest( fv );
    d_status->addStatusTest( absresid );
    d_status->addStatusTest( maxiters );
    d_status->addStatusTest( wrms );
    // Create nox parameter list
    d_nlParams             = Teuchos::rcp( new Teuchos::ParameterList );
    std::string solverType = nonlinear_db->getString( "solver" );
    if ( solverType == "JFNK" ) {
        d_nlParams->set( "Nonlinear Solver", "Line Search Based" );
    } else if ( solverType == "Anderson" ) {
        d_nlParams->set( "Nonlinear Solver", "Anderson Accelerated Fixed-Point" );
        int depth     = nonlinear_db->getIntegerWithDefault( "StorageDepth", 5 );
        double mixing = nonlinear_db->getDoubleWithDefault( "MixingParameter", 1.0 );
        d_nlParams->sublist( "Anderson Parameters" ).set( "Storage Depth", depth );
        d_nlParams->sublist( "Anderson Parameters" ).set( "Mixing Parameter", mixing );
        d_nlParams->sublist( "Anderson Parameters" )
            .sublist( "Preconditioning" )
            .set( "Precondition", d_precOp.get() != nullptr );
        Teuchos::RCP<NOX::StatusTest::RelativeNormF> relresid(
            new NOX::StatusTest::RelativeNormF( d_dMaxError ) );
        d_status->addStatusTest( relresid );
        Teuchos::RCP<AndersonStatusTest> andersonTest(
            new AMP::Solver::AndersonStatusTest( nonlinear_db ) );
        d_status->addStatusTest( andersonTest );
    }
    std::string lineSearchMethod =
        nonlinear_db->getStringWithDefault( "lineSearchMethod", "Polynomial" );
    d_nlParams->sublist( "Line Search" ).set( "Method", lineSearchMethod );
    d_nlParams->sublist( "Direction" )
        .sublist( "Newton" )
        .sublist( "Linear Solver" )
        .set( "Tolerance", linearRelativeTolerance );
    if ( params->d_prePostOperator.get() != nullptr ) {
        Teuchos::RefCountPtr<NOX::Abstract::PrePostOperator> prePostOperator(
            params->d_prePostOperator.get(),
            Teuchos::DeallocDelete<NOX::Abstract::PrePostOperator>(),
            false );
        d_nlParams->sublist( "Solver Options" )
            .set<Teuchos::RCP<NOX::Abstract::PrePostOperator>>( "User Defined Pre/Post Operator",
                                                                prePostOperator );
    }
    // Set the printing parameters in the "Printing" sublist
    Teuchos::ParameterList &printParams = d_nlParams->sublist( "Printing" );
    printParams.set( "Output Precision", 3 );
    printParams.set( "Output Processor", 0 );
    NOX::Utils::MsgType print_level = NOX::Utils::Error;
    if ( d_iDebugPrintInfoLevel >= 1 ) {
        print_level = static_cast<NOX::Utils::MsgType>(
            print_level + NOX::Utils::OuterIteration + NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::InnerIteration + NOX::Utils::Warning );
    } else if ( d_iDebugPrintInfoLevel >= 2 ) {
        print_level = static_cast<NOX::Utils::MsgType>(
            print_level + NOX::Utils::LinearSolverDetails + NOX::Utils::Parameters +
            NOX::Utils::Details + NOX::Utils::Debug + NOX::Utils::TestDetails + NOX::Utils::Error );
    }
    printParams.set( "Output Information", print_level );
}


/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
void TrilinosNOXSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                               AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    // PROFILE_START("solve");
    // Get thyra vectors
    auto initial = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( d_initialGuess ) );
    auto U = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::view( u ) );
    auto F = AMP::dynamic_pointer_cast<const AMP::LinearAlgebra::ThyraVector>(
        AMP::LinearAlgebra::ThyraVector::constView( f ) );
    NULL_USE( F );
    // Set the rhs for the thyra model
    d_thyraModel->setRhs( f );
    // Create the JFNK operator
    Teuchos::ParameterList printParams;
    auto jfnkParams = Teuchos::parameterList();
    jfnkParams->set( "Difference Type", "Forward" );
    jfnkParams->set( "Perturbation Algorithm", "KSP NOX 2001" );
    jfnkParams->set( "lambda", 1.0e-4 );
    Teuchos::RCP<NOX::Thyra::MatrixFreeJacobianOperator<double>> jfnkOp(
        new NOX::Thyra::MatrixFreeJacobianOperator<double>( printParams ) );
    jfnkOp->setParameterList( jfnkParams );
    if ( d_iDebugPrintInfoLevel >= 3 && d_comm.getRank() == 0 )
        jfnkParams->print( AMP::pout );
    // Create the NOX::Thyra::Group
    // Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group( initial->getVec(),
    // d_thyraModel ) );
    Teuchos::RCP<::Thyra::ModelEvaluator<double>> thyraModel(
        new NOX::MatrixFreeModelEvaluatorDecorator<double>( d_thyraModel ) );
    Teuchos::RCP<NOX::Thyra::Group> nox_group( new NOX::Thyra::Group(
        initial->getVec(), thyraModel, jfnkOp, d_lowsFactory, d_precOp, Teuchos::null ) );
    nox_group->setX( U->getVec() );
    nox_group->computeF();
    // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
    // This creates a circular dependency, so use a weak pointer.
    jfnkOp->setBaseEvaluationToNOXGroup( nox_group.create_weak() );
    // Create the solver
    d_solver = NOX::Solver::buildSolver( nox_group, d_status, d_nlParams );
    // Solve
    d_nlParams->print( AMP::pout );
    NOX::StatusTest::StatusType solvStatus = d_solver->solve();
    if ( solvStatus != NOX::StatusTest::Converged )
        AMP_ERROR( "Failed to solve" );
    // Copy the solution back to u
    const auto *tmp = dynamic_cast<const NOX::Thyra::Vector *>( &( nox_group->getX() ) );
    const auto *thyraVec =
        dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper *>( &( tmp->getThyraVector() ) );
    AMP_ASSERT( thyraVec != nullptr );
    AMP_ASSERT( thyraVec->numVecs() == 1 );
    u->copyVector( thyraVec->getVec( 0 ) );
    // PROFILE_STOP("solve");
}
} // namespace Solver
} // namespace AMP
