#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/vectors/ExternalVectorDeleter.h"
#include "AMP/vectors/petsc/ManagedPetscVector.h"
#include "ProfilerApp.h"

#include "petsc.h"
#include "petscksp.h"
#include "petscpc.h"

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


static inline PCSide getPCSide( const std::string &pc_side )
{
    PCSide PcSide = PC_RIGHT;
    if ( pc_side == "RIGHT" ) {
        PcSide = PC_RIGHT;
    } else if ( pc_side == "LEFT" ) {
        PcSide = PC_LEFT;
    } else if ( pc_side == "SYMMETRIC" ) {
        PcSide = PC_SYMMETRIC;
    } else {
        AMP_ERROR( "Unknown value for pc_type" );
    }
    return PcSide;
}


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
PetscKrylovSolver::PetscKrylovSolver()
{
    d_bKSPCreatedInternally = false;
    d_KrylovSolver          = nullptr;
}
PetscKrylovSolver::PetscKrylovSolver( AMP::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters )
{
    AMP_ASSERT( parameters.get() != nullptr );

    auto params = AMP::dynamic_pointer_cast<PetscKrylovSolverParameters>( parameters );
    AMP_ASSERT( params.get() != nullptr );

    // Create a default KrylovSolver
    d_bKSPCreatedInternally = true;
    KSPCreate( params->d_comm.getCommunicator(), &d_KrylovSolver );

    // Initialize
    initialize( parameters );
}


/****************************************************************
 *  De-constructor                                               *
 ****************************************************************/
PetscKrylovSolver::~PetscKrylovSolver()
{
    if ( d_bKSPCreatedInternally ) {
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        KSPDestroy( d_KrylovSolver );
#elif PETSC_VERSION_GE( 3, 2, 0 )
        KSPDestroy( &d_KrylovSolver );
#else
#error Not programmed for this version yet
#endif
    }
}


/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void PetscKrylovSolver::initialize( AMP::shared_ptr<SolverStrategyParameters> const params )
{
    auto parameters = AMP::dynamic_pointer_cast<PetscKrylovSolverParameters>( params );
    AMP_ASSERT( parameters.get() != nullptr );

    // the comm is set here of instead of the constructor because this routine
    // could be called by SNES directly
    d_comm = parameters->d_comm;
    AMP_ASSERT( !d_comm.isNull() );

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput( parameters->d_db );

    checkErr( KSPSetType( d_KrylovSolver, d_sKspType.c_str() ) );

    if ( d_KSPAppendOptionsPrefix != "" ) {
        KSPAppendOptionsPrefix( d_KrylovSolver, d_KSPAppendOptionsPrefix.c_str() );
        //      PCAppendOptionsPrefix(pc, d_KSPAppendOptionsPrefix.c_str());
    }

    if ( ( d_sKspType == "fgmres" ) || ( d_sKspType == "gmres" ) ) {
        checkErr( KSPGMRESSetRestart( d_KrylovSolver, d_iMaxKrylovDimension ) );
    }

    // Create the preconditioner
    PC pc;
    checkErr( KSPGetPC( d_KrylovSolver, &pc ) );
    if ( d_bUsesPreconditioner ) {
        if ( d_sPcType != "shell" ) {
            // the pointer to the preconditioner should be NULL if we are using a Petsc internal PC
            AMP_ASSERT( d_pPreconditioner.get() == nullptr );
            PCSetType( pc, d_sPcType.c_str() );
        } else {
            // for a shell preconditioner the user context is set to an instance of this class
            // and the setup and apply preconditioner functions for the PCSHELL
            // are set to static member functions of this class. By doing this we do not need to
            // introduce
            // static member functions into every SolverStrategy that might be used as a
            // preconditioner
            checkErr( PCSetType( pc, PCSHELL ) );
            checkErr( PCShellSetContext( pc, this ) );
            checkErr( PCShellSetSetUp( pc, PetscKrylovSolver::setupPreconditioner ) );
            checkErr( PCShellSetApply( pc, PetscKrylovSolver::applyPreconditioner ) );
        }
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        checkErr( KSPSetPreconditionerSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
#elif PETSC_VERSION_GE( 3, 2, 0 )
        checkErr( KSPSetPCSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
#else
#error Not programmed for this version yet
#endif
    } else {
        checkErr( PCSetType( pc, PCNONE ) );
    }

        // PetscTruth useZeroGuess = (d_bUseZeroInitialGuess) ? PETSC_TRUE : PETSC_FALSE;
        // ierr = KSPSetInitialGuessNonzero(d_KrylovSolver, useZeroGuess);

#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    PetscTruth useNonzeroGuess = ( !d_bUseZeroInitialGuess ) ? PETSC_TRUE : PETSC_FALSE;
#elif PETSC_VERSION_GE( 3, 2, 0 )
    PetscBool useNonzeroGuess = ( !d_bUseZeroInitialGuess ) ? PETSC_TRUE : PETSC_FALSE;
#else
#error Not programmed for this version yet
#endif
    checkErr( KSPSetInitialGuessNonzero( d_KrylovSolver, useNonzeroGuess ) );

    checkErr( KSPSetTolerances( d_KrylovSolver,
                                d_dRelativeTolerance,
                                d_dAbsoluteTolerance,
                                d_dDivergenceTolerance,
                                d_iMaxIterations ) );
    if ( d_bKSPCreatedInternally ) {
        checkErr( KSPSetFromOptions( d_KrylovSolver ) );
    }
    if ( d_PetscMonitor.get() != nullptr ) {
        // Add the monitor
        checkErr( KSPMonitorSet(
            d_KrylovSolver, PetscMonitor::monitorKSP, d_PetscMonitor.get(), PETSC_NULL ) );
    }
    // in this case we make the assumption we can access a PetscMat for now
    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }
}
// Function to get values from input
void PetscKrylovSolver::getFromInput( const AMP::shared_ptr<AMP::Database> &db )
{
    // fill this in
    std::string petscOptions = db->getStringWithDefault( "KSPOptions", "" );
    if ( petscOptions.find( "monitor" ) != std::string::npos ) {
        petscOptions = PetscMonitor::removeMonitor( petscOptions );
        d_PetscMonitor.reset( new PetscMonitor( d_comm ) );
    }
#if PETSC_VERSION_LT( 3, 3, 0 )
    PetscOptionsInsertString( petscOptions.c_str() );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
    PetscOptionsInsertString( PETSC_NULL, petscOptions.c_str() );
#else
    PetscOptions options;
    PetscOptionsInsertString( options, petscOptions.c_str() );
#error This does not seem right. The options are not used.  Check!!!
#endif


    d_sKspType             = db->getStringWithDefault( "ksp_type", "fgmres" );
    d_dRelativeTolerance   = db->getDoubleWithDefault( "relative_tolerance", 1.0e-9 );
    d_dAbsoluteTolerance   = db->getDoubleWithDefault( "absolute_tolerance", 1.0e-14 );
    d_dDivergenceTolerance = db->getDoubleWithDefault( "divergence_tolerance", 1.0e+03 );
    d_iMaxIterations       = db->getDoubleWithDefault( "max_iterations", 1000 );

    d_KSPAppendOptionsPrefix = db->getStringWithDefault( "KSPAppendOptionsPrefix", "" );

    if ( ( d_sKspType == "fgmres" ) || ( d_sKspType == "gmres" ) ) {
        d_iMaxKrylovDimension = db->getIntegerWithDefault( "max_krylov_dimension", 20 );
        d_sGmresOrthogonalizationAlgorithm =
            db->getStringWithDefault( "gmres_orthogonalization_algorithm", "modifiedgramschmidt" );
    }

    d_bUsesPreconditioner = db->getBoolWithDefault( "uses_preconditioner", false );

    if ( d_bUsesPreconditioner ) {
        if ( db->keyExists( "pc_type" ) ) {
            d_sPcType = db->getStringWithDefault( "pc_type", "none" );
        } else {
            // call error here
            AMP_ERROR( "pc_type does not exist" );
        }
        d_PcSide = db->getStringWithDefault( "pc_side", "RIGHT" );
    } else {
        d_sPcType = "none";
    }
}

/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
void PetscKrylovSolver::solve( AMP::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                               AMP::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    // fVecView and uVecView may be held in KSPSolve internals.
    // by declaring a temporary vector, we ensure that the KSPSolve
    // will be replaced by fVecView and uVecView before they are destroyed.
    AMP::LinearAlgebra::Vector::const_shared_ptr f_thisGetsAroundPETScSharedPtrIssue = fVecView;
    AMP::LinearAlgebra::Vector::shared_ptr u_thisGetsAroundPETScSharedPtrIssue       = uVecView;
#endif

// Get petsc views of the vectors
#if !( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    AMP::LinearAlgebra::Vector::const_shared_ptr fVecView;
    AMP::LinearAlgebra::Vector::shared_ptr uVecView;
#endif
    fVecView = AMP::LinearAlgebra::PetscVector::constView( f );
    uVecView = AMP::LinearAlgebra::PetscVector::view( u );

    // Check input vector states
    AMP_ASSERT(
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( fVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( fVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( uVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( uVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );

    if ( d_iDebugPrintInfoLevel > 1 ) {
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of rhs vector: " << f->L2Norm()
                  << std::endl;
    }
    Vec fVec = dynamic_pointer_cast<const AMP::LinearAlgebra::PetscVector>( fVecView )->getVec();
    Vec uVec = dynamic_pointer_cast<AMP::LinearAlgebra::PetscVector>( uVecView )->getVec();

    // Create the preconditioner and re-register the operator
    PC pc;
    checkErr( KSPGetPC( d_KrylovSolver, &pc ) );
    if ( d_bUsesPreconditioner ) {
        if ( d_sPcType != "shell" ) {
            // the pointer to the preconditioner should be NULL if we are using a Petsc internal PC
            AMP_ASSERT( d_pPreconditioner.get() == nullptr );
            PCSetType( pc, d_sPcType.c_str() );
        } else {
            // for a shell preconditioner the user context is set to an instance of this class
            // and the setup and apply preconditioner functions for the PCSHELL
            // are set to static member functions of this class. By doing this we do not need to
            // introduce
            // static member functions into every SolverStrategy that might be used as a
            // preconditioner
            checkErr( PCSetType( pc, PCSHELL ) );
            checkErr( PCShellSetContext( pc, this ) );
            checkErr( PCShellSetSetUp( pc, PetscKrylovSolver::setupPreconditioner ) );
            checkErr( PCShellSetApply( pc, PetscKrylovSolver::applyPreconditioner ) );
        }
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        checkErr( KSPSetPreconditionerSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
#elif PETSC_VERSION_GE( 3, 2, 0 )
        checkErr( KSPSetPCSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
#else
#error Not programmed for this version yet
#endif
    } else {
        checkErr( PCSetType( pc, PCNONE ) );
    }
    if ( d_pOperator.get() != nullptr ) {
        registerOperator( d_pOperator );
    }

    // This will replace any PETSc references to pointers we also track
    // After this, we are free to delet f_thisGetsAroundPETScSharedPtrIssue without memory leak.
    KSPSolve( d_KrylovSolver, fVec, uVec );

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution from KSP: " << u->L2Norm() << std::endl;
    }

// Reset the solvers
#if !( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
    KSPReset( d_KrylovSolver );
#endif
    PROFILE_STOP( "solve" );
}


/****************************************************************
 *  Function to set the KrylovSolver                             *
 ****************************************************************/
void PetscKrylovSolver::setKrylovSolver( KSP *ksp )
{
    if ( d_bKSPCreatedInternally ) {
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
        KSPDestroy( d_KrylovSolver );
#elif PETSC_VERSION_GE( 3, 2, 0 )
        KSPDestroy( &d_KrylovSolver );
#else
#error Not programmed for this version yet
#endif
    }
    d_bKSPCreatedInternally = false;
    d_KrylovSolver          = *ksp;
}


/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void PetscKrylovSolver::registerOperator( const AMP::shared_ptr<AMP::Operator::Operator> op )
{
    // in this case we make the assumption we can access a PetscMat for now
    AMP_ASSERT( op.get() != nullptr );

    d_pOperator = op;

    AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
        AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator.get() != nullptr );

    AMP::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>( linearOperator->getMatrix() );
    AMP_ASSERT( pMatrix.get() != nullptr );

    Mat mat;
    mat = pMatrix->getMat();

#if PETSC_VERSION_LE( 3, 2, 0 )
    KSPSetOperators( d_KrylovSolver, mat, mat, DIFFERENT_NONZERO_PATTERN );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
    KSPSetOperators( d_KrylovSolver, mat, mat );
#else
#error This version of PETSc is not supported.  Check!!!
#endif
}
void PetscKrylovSolver::resetOperator(
    const AMP::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator.get() != nullptr ) {
        d_pOperator->reset( params );
        AMP::shared_ptr<AMP::Operator::LinearOperator> linearOperator =
            AMP::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
        AMP_ASSERT( linearOperator.get() != nullptr );

        AMP::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix =
            AMP::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(
                linearOperator->getMatrix() );
        AMP_ASSERT( pMatrix.get() != nullptr );

        Mat mat;
        mat = pMatrix->getMat();


#if PETSC_VERSION_LE( 3, 2, 0 )
        KSPSetOperators( d_KrylovSolver, mat, mat, DIFFERENT_NONZERO_PATTERN );
#elif ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 7 )
        KSPSetOperators( d_KrylovSolver, mat, mat );
#else
#error This version of PETSc is not supported.  Check!!!
#endif // for now we will assume the pc is affected at every iteration
    }

    // should add a mechanism for the linear operator to provide updated parameters for the
    // preconditioner operator
    // though it's unclear where this might be necessary
    if ( d_pPreconditioner.get() != nullptr ) {
        d_pPreconditioner->resetOperator( params );
    }
}


/****************************************************************
 *  Function to setup the preconditioner                         *
 ****************************************************************/
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
int PetscKrylovSolver::setupPreconditioner( void * )
{
    int ierr = 0;
    // return( ((PetscKrylovSolver*)ctx)->getPreconditioner()->reset() );
    return ierr;
}
#elif PETSC_VERSION_GE( 3, 2, 0 )
PetscErrorCode PetscKrylovSolver::setupPreconditioner( PC pc )
{
    int ierr = 0;
    void *ctx = nullptr;
    ierr = PCShellGetContext( pc, &ctx );
    return ierr;
}
#else
#error Not programmed for this version yet
#endif


/****************************************************************
 *  Function to call the preconditioner                          *
 ****************************************************************/
#if ( PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0 )
PetscErrorCode PetscKrylovSolver::applyPreconditioner( void *ctx, Vec r, Vec z )
#elif PETSC_VERSION_GE( 3, 2, 0 )
PetscErrorCode PetscKrylovSolver::applyPreconditioner( PC pc, Vec r, Vec z )
#else
#error Not programmed for this version yet
#endif
{
    int ierr = 0;
#if PETSC_VERSION_GE( 3, 2, 0 )
    void *ctx;
    PCShellGetContext( pc, &ctx );
#endif
    AMP_ASSERT( ctx != nullptr );

    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_r(
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( r->data ),
        AMP::LinearAlgebra::ExternalVectorDeleter() );
    AMP::shared_ptr<AMP::LinearAlgebra::Vector> sp_z(
        reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>( z->data ),
        AMP::LinearAlgebra::ExternalVectorDeleter() );

    // Make sure the vectors are in a consistent state
    AMP_ASSERT(
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    AMP_ASSERT(
        ( sp_z->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::UNCHANGED ) ||
        ( sp_z->getUpdateStatus() == AMP::LinearAlgebra::Vector::UpdateState::LOCAL_CHANGED ) );
    sp_r->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
    sp_z->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );

    // these tests were helpful in finding a bug
    if ( ( (PetscKrylovSolver *) ctx )->getDebugPrintInfoLevel() > 5 ) {
        double norm = 0.0;
        VecNorm( r, NORM_2, &norm );
        double sp_r_norm = sp_r->L2Norm();
        AMP_ASSERT( AMP::Utilities::approx_equal( norm, sp_r_norm ) );
    }


    // Call the preconditioner
    AMP::shared_ptr<AMP::Solver::SolverStrategy> preconditioner =
        ( (PetscKrylovSolver *) ctx )->getPreconditioner();
    if ( preconditioner != nullptr ) {
        preconditioner->solve( sp_r, sp_z );
    } else {
        // Use the identity preconditioner
        sp_z->copyVector( sp_r );
    }

    // Check for nans (no communication necessary)
    double localNorm = sp_z->localL2Norm();
    AMP_INSIST( localNorm == localNorm, "NaNs detected in preconditioner" );

    // not sure why, but the state of sp_z is not updated
    // and petsc uses the cached norm
    auto firer = AMP::dynamic_pointer_cast<AMP::LinearAlgebra::DataChangeFirer>( sp_z );
    if ( firer )
        firer->fireDataChange();

    // these tests were helpful in finding a bug
    if ( ( (PetscKrylovSolver *) ctx )->getDebugPrintInfoLevel() > 5 ) {
        double norm = 0.0;
        AMP::pout << "L2 Norm of sp_z " << sp_z->L2Norm() << std::endl;
        VecNorm( z, NORM_2, &norm );
        AMP::pout << "L2 Norm of z " << norm << std::endl;
        AMP_ASSERT( norm == sp_z->L2Norm() );
    }

    return ( ierr );
}
}
} // namespace AMP
