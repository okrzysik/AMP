#include "PetscSNESSolver.h"

#include "ProfilerApp.h"
#include "utils/Utilities.h"

#include "matrices/petsc/PetscMatrix.h"
#include "operators/ColumnOperator.h"
#include "operators/LinearOperator.h"
#include "vectors/ExternalVectorDeleter.h"
#include "vectors/Vector.h"
#include "vectors/petsc/ManagedPetscVector.h"

#include "petscmat.h"
#include "petscsnes.h"


namespace AMP {
namespace Solver {


#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
static inline void checkErr( int ierr )
{
    AMP_INSIST( ierr == 0, "Petsc returned non-zero error code" );
}
#elif PETSC_VERSION_GE( 3, 2, 0 )
static inline void checkErr( PetscErrorCode ierr )
{
    AMP_INSIST( ierr == 0, "Petsc returned non-zero error code" );
}
#else
#error Not programmed for this version yet
#endif


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
PetscSNESSolver::PetscSNESSolver()
{
    d_bUsesJacobian                  = false;
    d_Jacobian                       = nullptr;
    d_SNESSolver                     = nullptr;
    d_sMFFDDifferencingStrategy      = MATMFFD_WP;
    d_dMFFDFunctionDifferencingError = PETSC_DEFAULT;
    d_pSolutionVector.reset();
    d_pResidualVector.reset();
    d_pKrylovSolver.reset();
}
PetscSNESSolver::PetscSNESSolver( AMP::shared_ptr<PetscSNESSolverParameters> parameters )
    : SolverStrategy( parameters )
{
    d_bUsesJacobian = false;
    d_Jacobian      = nullptr;
    d_SNESSolver    = nullptr;
    d_comm          = parameters->d_comm;
    d_pKrylovSolver = parameters->d_pKrylovSolver;
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
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    SNESDestroy( d_SNESSolver );
#elif PETSC_VERSION_GE( 3, 2, 0 )
    SNESDestroy( &d_SNESSolver );
#else
#error Not programmed for this version yet
#endif
    d_SNESSolver = nullptr;
}


/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void PetscSNESSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> params )
{
    PROFILE_START( "initialize" );

    auto parameters = AMP::dynamic_pointer_cast<PetscSNESSolverParameters>( params );
    getFromInput( parameters->d_db );

    // create the SNES solver
    if ( d_SNESSolver != nullptr ) {
        this->~PetscSNESSolver();
    }
    checkErr( SNESCreate( d_comm.getCommunicator(), &d_SNESSolver ) );

// set the type to line search, potentially modify this later to be from input
#if PETSC_VERSION_LE( 3, 2, 0 )
    checkErr( SNESSetType( d_SNESSolver, SNESLS ) );
#else
    checkErr( SNESSetType( d_SNESSolver, SNESNEWTONLS ) );
#endif
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
    if ( parameters->d_pInitialGuess.get() != nullptr ) {
        d_pSolutionVector = parameters->d_pInitialGuess;
    } else {
        AMP_INSIST( parameters->d_pInitialGuess.get() != nullptr,
                    "ERROR:: The initial guess has to "
                    "be provided through the "
                    "PetscSNESSolverParameters class" );
    }

    // if the krylov solver is initialized set the SNES pointer to it
    if ( d_pKrylovSolver.get() != nullptr ) {
        SNESSetKSP( d_SNESSolver, d_pKrylovSolver->getKrylovSolver() );
    } else {
        // initialize the Krylov solver correctly
        KSP kspSolver;
        // access the SNES internal pointer to KSP and get a pointer to KSP
        SNESGetKSP( d_SNESSolver, &kspSolver );

        const AMP::shared_ptr<AMP::Database> nonlinearSolverDb = parameters->d_db;

        if ( nonlinearSolverDb->keyExists( "LinearSolver" ) ) {
            d_pKrylovSolver.reset( new PetscKrylovSolver() );
            d_pKrylovSolver->setKrylovSolver( &kspSolver );
            PetscKrylovSolverParameters *krylovSolverParameters =
                new PetscKrylovSolverParameters( nonlinearSolverDb->getDatabase( "LinearSolver" ) );
            krylovSolverParameters->d_comm = d_comm;
            AMP::shared_ptr<SolverStrategyParameters> pKrylovSolverParameters(
                krylovSolverParameters );
            d_pKrylovSolver->initialize( pKrylovSolverParameters );
        } else {
            AMP_INSIST( d_pKrylovSolver.get() != nullptr,
                        "ERROR: The nonlinear solver database must "
                        "contain a database called LinearSolver" );
        }
    }

    if ( d_bEnableLineSearchPreCheck ) {
#if PETSC_VERSION_LE( 3, 2, 0 )
        checkErr(
            SNESLineSearchSetPreCheck( d_SNESSolver, PetscSNESSolver::lineSearchPreCheck, this ) );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
        SNESLineSearch snesLineSearch;
        SNESGetLineSearch( d_SNESSolver, &snesLineSearch );
        checkErr( SNESLineSearchSetPreCheck(
            snesLineSearch, PetscSNESSolver::lineSearchPreCheck, this ) );

#else
#error This version of PETSc is not supported.  Check!!!
#endif
    }

    checkErr( SNESSetFromOptions( d_SNESSolver ) );

    if ( d_PetscMonitor.get() != nullptr ) {
        // Add the monitor
        SNESMonitorSet( d_SNESSolver, PetscMonitor::monitorSNES, d_PetscMonitor.get(), PETSC_NULL );
    }
    PROFILE_STOP( "initialize" );
}
void PetscSNESSolver::getFromInput( const AMP::shared_ptr<AMP::Database> db )
{
    std::string petscOptions = db->getStringWithDefault( "SNESOptions", "" );
    d_PetscMonitor.reset();
    if ( petscOptions.find( "monitor" ) != std::string::npos ) {
        petscOptions = PetscMonitor::removeMonitor( petscOptions );
        d_PetscMonitor.reset( new PetscMonitor( d_comm ) );
    }
#if PETSC_VERSION_LE( 3, 2, 0 )
    PetscOptionsInsertString( petscOptions.c_str() );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
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
#else
#error This version of PETSc is not supported.  Check!!!
#endif

    d_bUsesJacobian = db->getBoolWithDefault( "usesJacobian", false );
    d_sMFFDDifferencingStrategy =
        db->getStringWithDefault( "MFFDDifferencingStrategy", MATMFFD_WP );
    d_dMFFDFunctionDifferencingError =
        db->getDoubleWithDefault( "MFFDFunctionDifferencingError", PETSC_DEFAULT );

    d_SNESAppendOptionsPrefix = db->getStringWithDefault( "SNESAppendOptionsPrefix", "" );

    if ( db->keyExists( "maximumFunctionEvals" ) )
        d_iMaximumFunctionEvals = db->getInteger( "maximumFunctionEvals" );

    if ( db->keyExists( "absoluteTolerance" ) )
        d_dAbsoluteTolerance = db->getDouble( "absoluteTolerance" );

    if ( db->keyExists( "relativeTolerance" ) )
        d_dRelativeTolerance = db->getDouble( "relativeTolerance" );

    if ( db->keyExists( "stepTolerance" ) )
        d_dStepTolerance = db->getDouble( "stepTolerance" );

    d_bEnableLineSearchPreCheck = db->getBoolWithDefault( "enableLineSearchPreCheck", false );

    if ( d_bEnableLineSearchPreCheck )
        d_iNumberOfLineSearchPreCheckAttempts =
            db->getIntegerWithDefault( "numberOfLineSearchPreCheckAttempts", 5 );

    d_bEnableMFFDBoundsCheck = db->getBoolWithDefault( "enableMFFDBoundsCheck", false );
    if ( d_bEnableMFFDBoundsCheck )
        d_operatorComponentToEnableBoundsCheck =
            db->getInteger( "operatorComponentToEnableBoundsCheck" );
}


/****************************************************************
 *  Apply                                                        *
 ****************************************************************/
PetscErrorCode PetscSNESSolver::apply( SNES, Vec x, Vec r, void *ctx )
{
    PROFILE_START( "apply" );
    int ierr = 0;

    auto *xvec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( x->data );
    auto *rvec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( r->data );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_x( xvec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_f;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_r( rvec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );

    if ( sp_f.get() != nullptr )
        sp_f->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    if ( sp_x.get() != nullptr )
        sp_x->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    AMP_ASSERT(
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    AMP::shared_ptr<AMP::Operator::Operator> op( ( (PetscSNESSolver *) ctx )->getOperator() );

    op->residual( sp_f, sp_x, sp_r );
    sp_r->scale( -1.0 );

    AMP_ASSERT(
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    PROFILE_STOP( "apply" );
    return ( ierr );
}


/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
void PetscSNESSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                             AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    if ( d_iDebugPrintInfoLevel > 2 )
        AMP::pout << "L2 Norm of u in PetscSNESSolver::solve before view " << u->L2Norm()
                  << std::endl;

    // Get petsc views of the vectors
    auto spRhs = AMP::LinearAlgebra::PetscVector::constView( f );
    auto spSol = AMP::LinearAlgebra::PetscVector::view( u );

// Create temporary copies of the petsc views
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    // This fixes a bug where a previous solve call creates and used views of a different vector,
    // which then are destroyed when the views of the new vectors are created, but petsc still
    // holds a copy of the original views until the new solve call is created
    d_refVectors.push_back( spRhs );
    d_refVectors.push_back( spSol );
    if ( d_pResidualVector != NULL )
        d_refVectors.push_back( AMP::LinearAlgebra::PetscVector::constView( d_pResidualVector ) );
#endif

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( spRhs->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( spRhs->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( spSol->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( spSol->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    if ( d_iDebugPrintInfoLevel > 2 )
        AMP::pout << "L2 Norm of u in PetscSNESSolver::solve after view " << spSol->L2Norm()
                  << std::endl;

    // if the dynamic cast yielded a valid pointer
    if ( spSol.get() != nullptr ) {

        Vec x = dynamic_pointer_cast<const AMP::LinearAlgebra::PetscVector>( spSol )->getVec();

        Vec b = PETSC_NULL;
        if ( spRhs.get() != nullptr ) {
            b = dynamic_pointer_cast<const AMP::LinearAlgebra::PetscVector>( spRhs )->getVec();
            setSNESFunction( spRhs );
        }

        // Set the jacobian
        if ( !d_bUsesJacobian ) {
            if ( d_Jacobian != nullptr ) {
                PETSC::matDestroy( &d_Jacobian );
                d_Jacobian = nullptr;
            }
            checkErr( MatCreateSNESMF( d_SNESSolver, &d_Jacobian ) );
            checkErr(
                MatMFFDSetType( d_Jacobian, (MatMFFDType) d_sMFFDDifferencingStrategy.c_str() ) );
            checkErr( MatMFFDSetFunctionError( d_Jacobian, d_dMFFDFunctionDifferencingError ) );
            if ( d_bEnableMFFDBoundsCheck ) {
                checkErr( MatMFFDSetCheckh( d_Jacobian, PetscSNESSolver::mffdCheckBounds, this ) );
            }
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
            checkErr( MatMFFDSetFromOptions( d_Jacobian ) );
#elif PETSC_VERSION_GE( 3, 2, 0 )
            checkErr( MatSetFromOptions( d_Jacobian ) );
#else
#error Not programmed for this version yet
#endif
        } else {
            AMP::shared_ptr<AMP::Operator::LinearOperator> linearOp =
                AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>(
                    d_pKrylovSolver->getOperator() );
            if ( linearOp.get() != nullptr ) {
                AMP::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix =
                    AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(
                        linearOp->getMatrix() );
                AMP_ASSERT( pMatrix.get() != nullptr );
                d_Jacobian = pMatrix->getMat();
            } else {
                AMP_INSIST( linearOp.get() != nullptr,
                            "ERROR: The LinearOperator pointer in the PetscKrylovSolver is NULL" );
            }
        }
        AMP::shared_ptr<AMP::Solver::SolverStrategy> pcSolver =
            d_pKrylovSolver->getPreconditioner();
        Mat PCJacobian = d_Jacobian;
        if ( pcSolver.get() != nullptr ) {
            AMP::shared_ptr<AMP::Operator::LinearOperator> linearOp =
                AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( pcSolver->getOperator() );
            if ( linearOp.get() != nullptr ) {
                AMP::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix =
                    AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(
                        linearOp->getMatrix() );
                if ( pMatrix.get() != nullptr ) {
                    PCJacobian = pMatrix->getMat();
                }
            }
        }

        checkErr( SNESSetJacobian(
            d_SNESSolver, d_Jacobian, PCJacobian, PetscSNESSolver::setJacobian, this ) );

        // Solve
        PROFILE_START( "petsc-SNESSolve" );
        checkErr( SNESSolve( d_SNESSolver, b, x ) );
        PROFILE_STOP( "petsc-SNESSolve" );
    } else {
        AMP_INSIST( spSol.get() != nullptr,
                    "ERROR: Currently the SNES Solver can only be used with a Petsc_Vector, the "
                    "supplied Vector does not appear to belong to this class" );
    }

// Reset the solvers
#if !( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    SNESReset( d_SNESSolver );
#endif

// Delete any copies of the reference vectors that we can
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    std::list<AMP::LinearAlgebra::Vector::const_shared_ptr>::iterator it = d_refVectors.begin();
    while ( it != d_refVectors.end() ) {
        const AMP::LinearAlgebra::PetscVector *view =
            dynamic_cast<const AMP::LinearAlgebra::PetscVector *>( it->get() );
        AMP_ASSERT( view != NULL );
        if ( view->petscHoldsView() )
            ++it;
        else
            it = d_refVectors.erase( it );
    }
#endif

    spRhs.reset();
    spSol.reset();

    u->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    PROFILE_STOP( "solve" );
}


/****************************************************************
 *  setJacobian                                                  *
 ****************************************************************/
#if PETSC_VERSION_LE( 3, 2, 0 )
PetscErrorCode PetscSNESSolver::setJacobian( SNES, Vec x, Mat *A, Mat *, MatStructure *, void *ctx )
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
PetscErrorCode PetscSNESSolver::setJacobian( SNES, Vec x, Mat A, Mat, void *ctx )
#else
#error This version of PETSc is not supported.  Check!!!
#endif
{
    PROFILE_START( "setJacobian" );
    int ierr           = 0;
    auto *pSNESSolver  = (PetscSNESSolver *) ctx;
    bool bUsesJacobian = pSNESSolver->getUsesJacobian();

    if ( !bUsesJacobian ) {
#if PETSC_VERSION_LE( 3, 2, 0 )
        ierr = MatAssemblyBegin( *A, MAT_FINAL_ASSEMBLY );
        ierr = MatAssemblyEnd( *A, MAT_FINAL_ASSEMBLY );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
        ierr = MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
        ierr = MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
#else
#error This version of PETSc is not supported.  Check!!!
#endif
    }

    auto *pVecShell = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( x->data );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> pSolution(
        pVecShell, AMP::LinearAlgebra::ExternalVectorDeleter() );

    AMP::shared_ptr<AMP::Operator::Operator> op = pSNESSolver->getOperator();
    AMP::shared_ptr<AMP::Operator::OperatorParameters> op_parameters =
        op->getParameters( "Jacobian", pSolution );
    AMP::shared_ptr<PetscKrylovSolver> pKrylovSolver = pSNESSolver->getKrylovSolver();
    pKrylovSolver->resetOperator( op_parameters );

    PROFILE_STOP( "setJacobian" );
    return ierr;
}


/****************************************************************
 *  Check if the vector is valid                                 *
 ****************************************************************/
bool PetscSNESSolver::isVectorValid( AMP::shared_ptr<AMP::Operator::Operator> &op,
                                     AMP::LinearAlgebra::Vector::shared_ptr &v,
                                     AMP_MPI comm )
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
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
PetscErrorCode
PetscSNESSolver::lineSearchPreCheck( SNES, Vec x, Vec y, void *checkctx, PetscTruth *changed_y )
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 2 )
PetscErrorCode
PetscSNESSolver::lineSearchPreCheck( SNES, Vec x, Vec y, void *checkctx, PetscBool *changed_y )
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
PetscErrorCode PetscSNESSolver::lineSearchPreCheck(
    SNESLineSearch, Vec x, Vec y, PetscBool *changed_y, void *checkctx )
#else
#error Not programmed for this version yet
#endif
{
    int ierr          = 1;
    auto *pSNESSolver = (PetscSNESSolver *) checkctx;

    AMP::shared_ptr<AMP::Operator::Operator> pOperator         = pSNESSolver->getOperator();
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> pScratchVector = pSNESSolver->getScratchVector();

    auto *xvec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( x->data );
    auto *yvec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( y->data );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_x( xvec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_y( yvec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );

    pScratchVector->add( sp_x, sp_y );

    if ( isVectorValid( pOperator, pScratchVector, xvec->getComm() ) ) {
        *changed_y = PETSC_FALSE;
        ierr       = 0;
    } else {
        int iNumberOfLineSearchPreCheckAttempts =
            pSNESSolver->getNumberOfLineSearchPreCheckAttempts();
        AMP::shared_ptr<AMP::Operator::ColumnOperator> pColumnOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( pOperator );
        if ( pColumnOperator.get() != nullptr ) {
            for ( int i = 0; i < iNumberOfLineSearchPreCheckAttempts; i++ ) {
                AMP::pout << "Attempting to scale search, attempt number " << i << std::endl;
                double lambda = 0.5;
                sp_y->scale( lambda, sp_y );
                pScratchVector->add( sp_x, sp_y );
                if ( isVectorValid( pOperator, pScratchVector, xvec->getComm() ) ) {
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
    auto *pSNESSolver                                      = (PetscSNESSolver *) checkctx;
    AMP::shared_ptr<AMP::Operator::Operator> pSNESOperator = pSNESSolver->getOperator();
    AMP::shared_ptr<AMP::Operator::Operator> pOperator;
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> pScratchVector = pSNESSolver->getScratchVector();

    auto *uvec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( U->data );
    auto *avec = reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( a->data );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_u( uvec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_a( avec,
                                                      AMP::LinearAlgebra::ExternalVectorDeleter() );

    // check for column operators
    AMP::shared_ptr<AMP::Operator::ColumnOperator> pColumnOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::ColumnOperator>( pSNESOperator );
    if ( pColumnOperator.get() != nullptr ) {
        pOperator = pColumnOperator->getOperator( pSNESSolver->getBoundsCheckComponent() );
    } else {
        pOperator = pSNESOperator;
    }

    AMP::LinearAlgebra::Variable::shared_ptr opVar = pOperator->getOutputVariable();
    AMP::LinearAlgebra::Vector::shared_ptr scv = pScratchVector->subsetVectorForVariable( opVar );
    AMP::LinearAlgebra::Vector::shared_ptr uv  = sp_u->subsetVectorForVariable( opVar );
    AMP::LinearAlgebra::Vector::shared_ptr av  = sp_a->subsetVectorForVariable( opVar );

    scv->axpy( *h, av, uv );

    // the code below is only valid for ensuring positivity
    // will do for now
    if ( isVectorValid( pOperator, scv, uvec->getComm() ) ) {
        double minVal = PetscAbsScalar( ( *h ) * 1.01 );
        scv->divide( uv, av );
        scv->abs( scv );
        minVal = std::min( scv->min(), minVal );
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


void PetscSNESSolver::setSNESFunction( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> rhs )
{
    AMP_INSIST( rhs.get() != nullptr,
                "ERROR: PetscSNESSolver::setSNESFunction needs a non NULL rhs vector argument" );

    // Create new residual and scratch vectors
    d_pResidualVector = rhs->cloneVector();
    d_pScratchVector  = d_pResidualVector->cloneVector();

    // set the function evaluation routine to a static member of this class which acts as a wrapper
    AMP::LinearAlgebra::Vector::shared_ptr petscVec =
        AMP::LinearAlgebra::PetscVector::view( d_pResidualVector );
    AMP_INSIST( petscVec.get() != nullptr,
                "ERROR: Currently the SNES Solver can only be used with a Petsc_Vector, the "
                "supplied Vector does not appear to belong to this class" );
    Vec residualVector =
        dynamic_pointer_cast<AMP::LinearAlgebra::PetscVector>( petscVec )->getVec();
    SNESSetFunction( d_SNESSolver, residualVector, PetscSNESSolver::apply, (void *) this );
}

void PetscSNESSolver::setInitialGuess( AMP::shared_ptr<AMP::LinearAlgebra::Vector> initialGuess )
{
    d_pSolutionVector->copyVector( initialGuess );
}


} // Solver
} // namespace AMP
