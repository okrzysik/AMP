#include "AMP/solvers/petsc/PetscSNESSolver.h"
#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/operators/ColumnOperator.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Vector.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "ProfilerApp.h"

#include "petsc/private/vecimpl.h"
#include "petscmat.h"
#include "petscsnes.h"


namespace AMP {
namespace Solver {


#if PETSC_VERSION_LT( 3, 7, 5 )
    #error AMP only supports PETSc 3.7.5 or greater
#endif


static inline void checkErr( PetscErrorCode ierr )
{
    AMP_INSIST( ierr == 0, "Petsc returned non-zero error code" );
}


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
PetscSNESSolver::PetscSNESSolver()
    : d_bUsesJacobian( false ),
      d_bEnableLineSearchPreCheck( false ),
      d_bEnableMFFDBoundsCheck( false ),
      d_iMaximumFunctionEvals( 0 ),
      d_iNumberOfLineSearchPreCheckAttempts( 0 ),
      d_operatorComponentToEnableBoundsCheck( 0 ),
      d_dStepTolerance( 0 ),
      d_sMFFDDifferencingStrategy( MATMFFD_WP ),
      d_dMFFDFunctionDifferencingError( PETSC_DEFAULT ),
      d_comm( AMP_COMM_NULL ),
      d_pSolutionVector( nullptr ),
      d_pResidualVector( nullptr ),
      d_SNESSolver( nullptr ),
      d_Jacobian( nullptr ),
      d_pKrylovSolver( nullptr )
{
}
PetscSNESSolver::PetscSNESSolver( std::shared_ptr<PetscSNESSolverParameters> parameters )
    : SolverStrategy( parameters ),
      d_bUsesJacobian( false ),
      d_bEnableLineSearchPreCheck( false ),
      d_bEnableMFFDBoundsCheck( false ),
      d_iMaximumFunctionEvals( 0 ),
      d_iNumberOfLineSearchPreCheckAttempts( 0 ),
      d_operatorComponentToEnableBoundsCheck( 0 ),
      d_dStepTolerance( 0 ),
      d_sMFFDDifferencingStrategy( MATMFFD_WP ),
      d_dMFFDFunctionDifferencingError( PETSC_DEFAULT ),
      d_comm( parameters->d_comm ),
      d_pSolutionVector( nullptr ),
      d_pResidualVector( nullptr ),
      d_SNESSolver( nullptr ),
      d_Jacobian( nullptr ),
      d_pKrylovSolver( parameters->d_pKrylovSolver )
{
    initialize( parameters );
}


/****************************************************************
 *  De-constructor                                               *
 ****************************************************************/
PetscSNESSolver::~PetscSNESSolver()
{
    // when we are using Matrix free delete the MF PETSc Jacobian
    if ( ( !d_bUsesJacobian ) && ( d_Jacobian != nullptr ) ) {
        PETSC::matDestroy( &d_Jacobian );
        d_Jacobian = nullptr;
    }
    SNESMonitorCancel( d_SNESSolver );
    SNESDestroy( &d_SNESSolver );
    d_SNESSolver = nullptr;
}


/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void PetscSNESSolver::initialize( std::shared_ptr<const SolverStrategyParameters> params )
{
    PROFILE_START( "initialize" );

    auto parameters = std::dynamic_pointer_cast<const PetscSNESSolverParameters>( params );
    getFromInput( parameters->d_db );

    // create the SNES solver
    if ( d_SNESSolver != nullptr ) {
        this->~PetscSNESSolver();
    }
    checkErr( SNESCreate( d_comm.getCommunicator(), &d_SNESSolver ) );

    // set the type to line search, potentially modify this later to be from input
    checkErr( SNESSetType( d_SNESSolver, SNESNEWTONLS ) );
    checkErr( SNESSetApplicationContext( d_SNESSolver, this ) );
    checkErr( SNESSetTolerances( d_SNESSolver,
                                 d_dAbsoluteTolerance,
                                 d_dRelativeTolerance,
                                 d_dStepTolerance,
                                 d_iMaxIterations,
                                 d_iMaximumFunctionEvals ) );
    if ( d_SNESAppendOptionsPrefix != "" )
        SNESAppendOptionsPrefix( d_SNESSolver, d_SNESAppendOptionsPrefix.c_str() );

    // if the initial guess is non-zero set the vectors accordingly
    if ( parameters->d_pInitialGuess ) {
        d_pSolutionVector = parameters->d_pInitialGuess;
    } else {
        AMP_INSIST( parameters->d_pInitialGuess,
                    "ERROR:: The initial guess has to "
                    "be provided through the "
                    "PetscSNESSolverParameters class" );
    }

    // if the krylov solver is initialized set the SNES pointer to it
    if ( d_pKrylovSolver ) {
        SNESSetKSP( d_SNESSolver, d_pKrylovSolver->getKrylovSolver() );
    } else {
        // initialize the Krylov solver correctly
        KSP kspSolver;
        // access the SNES internal pointer to KSP and get a pointer to KSP
        SNESGetKSP( d_SNESSolver, &kspSolver );

        auto nonlinearSolverDb = parameters->d_db;

        if ( nonlinearSolverDb->keyExists( "LinearSolver" ) ) {
            d_pKrylovSolver.reset( new PetscKrylovSolver() );
            d_pKrylovSolver->setKrylovSolver( &kspSolver );
            auto params2 = std::make_shared<PetscKrylovSolverParameters>(
                nonlinearSolverDb->getDatabase( "LinearSolver" ) );
            params2->d_comm = d_comm;
            d_pKrylovSolver->initialize( params2 );
        } else {
            AMP_INSIST( d_pKrylovSolver,
                        "ERROR: The nonlinear solver database must "
                        "contain a database called LinearSolver" );
        }
    }

    if ( d_bEnableLineSearchPreCheck ) {
        SNESLineSearch snesLineSearch;
        SNESGetLineSearch( d_SNESSolver, &snesLineSearch );
        checkErr( SNESLineSearchSetPreCheck(
            snesLineSearch, PetscSNESSolver::lineSearchPreCheck, this ) );
    }

    checkErr( SNESSetFromOptions( d_SNESSolver ) );

    if ( d_PetscMonitor ) {
        // Add the monitor
        SNESMonitorSet( d_SNESSolver, PetscMonitor::monitorSNES, d_PetscMonitor.get(), PETSC_NULL );
    }
    PROFILE_STOP( "initialize" );
}
void PetscSNESSolver::getFromInput( std::shared_ptr<const AMP::Database> db )
{
    std::string petscOptions = db->getWithDefault<std::string>( "SNESOptions", "" );
    d_PetscMonitor.reset();
    if ( petscOptions.find( "monitor" ) != std::string::npos ) {
        petscOptions = PetscMonitor::removeMonitor( petscOptions );
        d_PetscMonitor.reset( new PetscMonitor( d_comm ) );
    }
    // if the solver type is specified as 'ls' replace by 'newtonls'
    // this hack is brittle and will easily fail if the string is
    // not matched exactly in the input. Temporary fix for now
    std::string solverTypeStr( "-snes_type ls" );
    auto pos = petscOptions.find( solverTypeStr );
    if ( pos != std::string::npos ) {
        petscOptions.erase( pos, solverTypeStr.length() );
        petscOptions += " -snes_type newtonls";
    }

    PetscOptionsInsertString( PETSC_NULL, petscOptions.c_str() );

    d_bUsesJacobian = db->getWithDefault<bool>( "usesJacobian", false );
    d_sMFFDDifferencingStrategy =
        db->getWithDefault<std::string>( "MFFDDifferencingStrategy", MATMFFD_WP );
    d_dMFFDFunctionDifferencingError =
        db->getWithDefault<double>( "MFFDFunctionDifferencingError", PETSC_DEFAULT );

    d_SNESAppendOptionsPrefix = db->getWithDefault<std::string>( "SNESAppendOptionsPrefix", "" );

    if ( db->keyExists( "maximumFunctionEvals" ) )
        d_iMaximumFunctionEvals = db->getScalar<int>( "maximumFunctionEvals" );

    if ( db->keyExists( "stepTolerance" ) )
        d_dStepTolerance = db->getScalar<double>( "stepTolerance" );

    d_bEnableLineSearchPreCheck = db->getWithDefault<bool>( "enableLineSearchPreCheck", false );

    if ( d_bEnableLineSearchPreCheck )
        d_iNumberOfLineSearchPreCheckAttempts =
            db->getWithDefault<int>( "numberOfLineSearchPreCheckAttempts", 5 );

    d_bEnableMFFDBoundsCheck = db->getWithDefault<bool>( "enableMFFDBoundsCheck", false );
    if ( d_bEnableMFFDBoundsCheck )
        d_operatorComponentToEnableBoundsCheck =
            db->getScalar<int>( "operatorComponentToEnableBoundsCheck" );
}


/****************************************************************
 *  Apply                                                        *
 ****************************************************************/
PetscErrorCode PetscSNESSolver::apply( SNES, Vec x, Vec r, void *ctx )
{
    PROFILE_START( "apply" );
    int ierr = 0;

    auto sp_x = PETSC::getAMP( x );
    auto sp_r = PETSC::getAMP( r );

    std::shared_ptr<AMP::LinearAlgebra::Vector> sp_f;
    if ( sp_f )
        sp_f->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    if ( sp_x )
        sp_x->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    sp_r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    auto *pSNESSolver = reinterpret_cast<PetscSNESSolver *>( ctx );
    std::shared_ptr<AMP::Operator::Operator> op( pSNESSolver->getOperator() );

    op->residual( sp_f, sp_x, sp_r );
    sp_r->scale( -1.0 );
    sp_r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    PROFILE_STOP( "apply" );
    return ( ierr );
}


/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
void PetscSNESSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                             std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    if ( d_iDebugPrintInfoLevel > 2 )
        AMP::pout << "L2 Norm of u in PetscSNESSolver::solve before view " << u->L2Norm()
                  << std::endl;

    // Get petsc views of the vectors
    auto spRhs = AMP::LinearAlgebra::PetscVector::constView( f );
    auto spSol = AMP::LinearAlgebra::PetscVector::view( u );
    AMP_ASSERT( spSol );

    // Check input vector states
    using UpdateState = AMP::LinearAlgebra::VectorData::UpdateState;
    AMP_ASSERT( ( f->getUpdateStatus() == UpdateState::UNCHANGED ) ||
                ( f->getUpdateStatus() == UpdateState::LOCAL_CHANGED ) );
    u->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    if ( d_iDebugPrintInfoLevel > 2 )
        AMP::pout << "L2 Norm of u in PetscSNESSolver::solve after view " << u->L2Norm()
                  << std::endl;

    Vec x = spSol->getVec();

    Vec b = PETSC_NULL;
    if ( spRhs ) {
        b = spRhs->getVec();
        setSNESFunction( spRhs->getManagedVec() );
    }

    // Set the jacobian
    std::shared_ptr<AMP::LinearAlgebra::PetscMatrix> view1;
    if ( !d_bUsesJacobian ) {
        if ( d_Jacobian ) {
            PETSC::matDestroy( &d_Jacobian );
            d_Jacobian = nullptr;
        }
        checkErr( MatCreateSNESMF( d_SNESSolver, &d_Jacobian ) );
        checkErr( MatMFFDSetType( d_Jacobian, (MatMFFDType) d_sMFFDDifferencingStrategy.c_str() ) );
        checkErr( MatMFFDSetFunctionError( d_Jacobian, d_dMFFDFunctionDifferencingError ) );
        if ( d_bEnableMFFDBoundsCheck ) {
            checkErr( MatMFFDSetCheckh( d_Jacobian, PetscSNESSolver::mffdCheckBounds, this ) );
        }
        checkErr( MatSetFromOptions( d_Jacobian ) );
    } else {
        auto linearOp = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>(
            d_pKrylovSolver->getOperator() );
        if ( linearOp ) {
            view1      = AMP::LinearAlgebra::PetscMatrix::view( linearOp->getMatrix() );
            d_Jacobian = view1->getMat();
        } else {
            AMP_INSIST( linearOp,
                        "ERROR: The LinearOperator pointer in the PetscKrylovSolver is NULL" );
        }
    }
    auto pcSolver  = d_pKrylovSolver->getPreconditioner();
    Mat PCJacobian = d_Jacobian;
    std::shared_ptr<AMP::LinearAlgebra::PetscMatrix> view2;
    if ( pcSolver ) {
        auto linearOp =
            std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( pcSolver->getOperator() );
        if ( linearOp ) {
            auto matrix = linearOp->getMatrix();
            if ( matrix ) {
                view2      = AMP::LinearAlgebra::PetscMatrix::view( matrix );
                PCJacobian = view2->getMat();
            }
        }
    }

    checkErr( SNESSetJacobian(
        d_SNESSolver, d_Jacobian, PCJacobian, PetscSNESSolver::setJacobian, this ) );

    // Solve
    PROFILE_START( "petsc-SNESSolve" );
    checkErr( SNESSolve( d_SNESSolver, b, x ) );
    PROFILE_STOP( "petsc-SNESSolve" );

    // Reset the solvers
    SNESReset( d_SNESSolver );

    spRhs.reset();
    spSol.reset();

    u->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    PROFILE_STOP( "solve" );
}


/****************************************************************
 *  setJacobian                                                  *
 ****************************************************************/
PetscErrorCode PetscSNESSolver::setJacobian( SNES, Vec x, Mat A, Mat, void *ctx )
{
    PROFILE_START( "setJacobian" );
    int ierr           = 0;
    auto *pSNESSolver  = reinterpret_cast<PetscSNESSolver *>( ctx );
    bool bUsesJacobian = pSNESSolver->getUsesJacobian();

    if ( !bUsesJacobian ) {
        ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
        ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
    }

    auto pSolution     = PETSC::getAMP( x );
    auto op            = pSNESSolver->getOperator();
    auto op_parameters = op->getParameters( "Jacobian", pSolution );
    auto pKrylovSolver = pSNESSolver->getKrylovSolver();
    pKrylovSolver->resetOperator( op_parameters );

    PROFILE_STOP( "setJacobian" );
    return ierr;
}


/****************************************************************
 *  Check if the vector is valid                                 *
 ****************************************************************/
bool PetscSNESSolver::isVectorValid( std::shared_ptr<AMP::Operator::Operator> &op,
                                     AMP::LinearAlgebra::Vector::shared_ptr &v,
                                     const AMP_MPI &comm )
{
    bool retVal = false;
    int msg     = op->isValidInput( v ) ? 1 : 0;
    int result  = comm.minReduce( msg );
    retVal      = ( result == 1 );
    return retVal;
}


/****************************************************************
 *  Linesearch precheck                                          *
 ****************************************************************/
PetscErrorCode PetscSNESSolver::lineSearchPreCheck(
    SNESLineSearch, Vec x, Vec y, PetscBool *changed_y, void *checkctx )
{
    int ierr          = 1;
    auto *pSNESSolver = reinterpret_cast<PetscSNESSolver *>( checkctx );

    auto pOperator      = pSNESSolver->getOperator();
    auto pScratchVector = pSNESSolver->getScratchVector();

    auto sp_x = PETSC::getAMP( x );
    auto sp_y = PETSC::getAMP( y );

    pScratchVector->add( *sp_x, *sp_y );

    if ( isVectorValid( pOperator, pScratchVector, sp_x->getComm() ) ) {
        *changed_y = PETSC_FALSE;
        ierr       = 0;
    } else {
        int N_line = pSNESSolver->getNumberOfLineSearchPreCheckAttempts();
        auto pColumnOperator =
            std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( pOperator );
        if ( pColumnOperator ) {
            for ( int i = 0; i < N_line; i++ ) {
                AMP::pout << "Attempting to scale search, attempt number " << i << std::endl;
                double lambda = 0.5;
                sp_y->scale( lambda, *sp_y );
                pScratchVector->add( *sp_x, *sp_y );
                if ( isVectorValid( pOperator, pScratchVector, sp_x->getComm() ) ) {
                    ierr       = 0;
                    *changed_y = PETSC_TRUE;
                    break;
                } else {
                    lambda = lambda / 2.0;
                }
            }
        }
    }
    return ierr;
}


PetscErrorCode PetscSNESSolver::mffdCheckBounds( void *checkctx, Vec U, Vec a, PetscScalar *h )
{
    auto *pSNESSolver  = reinterpret_cast<PetscSNESSolver *>( checkctx );
    auto pSNESOperator = pSNESSolver->getOperator();
    std::shared_ptr<AMP::Operator::Operator> pOperator;
    auto pScratchVector = pSNESSolver->getScratchVector();

    auto sp_u = PETSC::getAMP( U );
    auto sp_a = PETSC::getAMP( a );

    // check for column operators
    auto pColumnOperator =
        std::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( pSNESOperator );
    if ( pColumnOperator ) {
        pOperator = pColumnOperator->getOperator( pSNESSolver->getBoundsCheckComponent() );
    } else {
        pOperator = pSNESOperator;
    }

    auto opVar = pOperator->getOutputVariable();
    auto scv   = pScratchVector->subsetVectorForVariable( opVar );
    auto uv    = sp_u->subsetVectorForVariable( opVar );
    auto av    = sp_a->subsetVectorForVariable( opVar );

    scv->axpy( *h, *av, *uv );

    // the code below is only valid for ensuring positivity
    // will do for now
    if ( isVectorValid( pOperator, scv, sp_u->getComm() ) ) {
        double minVal = PetscAbsScalar( ( *h ) * 1.01 );
        scv->divide( *uv, *av );
        scv->abs( *scv );
        minVal = std::min( static_cast<double>( scv->min() ), minVal );
        if ( minVal <= PetscAbsScalar( *h ) ) {
            AMP::pout << "Scaling h back from  " << ( *h ) << " to " << 0.99 * minVal << std::endl;
            if ( PetscRealPart( *h ) > 0.0 )
                *h = 0.99 * minVal;
            else
                *h = -0.99 * minVal;
        }
    }

    return ( 0 );
}


void PetscSNESSolver::setSNESFunction( std::shared_ptr<const AMP::LinearAlgebra::Vector> rhs )
{
    AMP_INSIST( rhs,
                "ERROR: PetscSNESSolver::setSNESFunction needs a non NULL rhs vector argument" );

    // Create new residual and scratch vectors
    d_pResidualVector = rhs->cloneVector();
    d_pScratchVector  = d_pResidualVector->cloneVector();

    // set the function evaluation routine to a static member of this class which acts as a wrapper
    auto petscVec = AMP::LinearAlgebra::PetscVector::view( d_pResidualVector );
    AMP_INSIST( petscVec,
                "ERROR: Currently the SNES Solver can only be used with a Petsc_Vector, "
                "the supplied Vector does not appear to belong to this class" );
    Vec residualVector = petscVec->getVec();
    SNESSetFunction( d_SNESSolver, residualVector, PetscSNESSolver::apply, (void *) this );
}

void PetscSNESSolver::setInitialGuess( std::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess )
{
    d_pSolutionVector->copyVector( initialGuess );
}


} // namespace Solver
} // namespace AMP
