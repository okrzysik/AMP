#include "AMP/solvers/petsc/PetscKrylovSolver.h"
#include "AMP/matrices/Matrix.h"
#include "AMP/matrices/petsc/PetscMatrix.h"
#include "AMP/operators/LinearOperator.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/petsc/PetscHelpers.h"
#include "AMP/vectors/petsc/PetscVector.h"

#include "ProfilerApp.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscksp.h"
#include "petscpc.h"

namespace AMP {
namespace Solver {


#if PETSC_VERSION_LT( 3, 7, 5 )
#error AMP only supports PETSc 3.7.5 or greater
#endif


static inline void checkErr( PetscErrorCode ierr )
{
    AMP_INSIST( ierr == 0, "Petsc returned non-zero error code" );
}


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


static inline void checkUpdateStatus( std::shared_ptr<const AMP::LinearAlgebra::Vector> x )
{
    auto status = x->getUpdateStatus();
    AMP_ASSERT( ( status == AMP::LinearAlgebra::VectorData::UpdateState::UNCHANGED ) ||
                ( status == AMP::LinearAlgebra::VectorData::UpdateState::LOCAL_CHANGED ) );
}


/****************************************************************
 *  Constructors                                                 *
 ****************************************************************/
PetscKrylovSolver::PetscKrylovSolver()
    : d_dDivergenceTolerance( 0 ),
      d_bKSPCreatedInternally( false ),
      d_bUsesPreconditioner( false ),
      d_iMaxKrylovDimension( 0 ),
      d_KrylovSolver( nullptr )
{
}
PetscKrylovSolver::PetscKrylovSolver( std::shared_ptr<SolverStrategyParameters> parameters )
    : SolverStrategy( parameters ),
      d_dDivergenceTolerance( 0 ),
      d_bKSPCreatedInternally( false ),
      d_bUsesPreconditioner( false ),
      d_iMaxKrylovDimension( 0 ),
      d_KrylovSolver( nullptr )
{
    AMP_ASSERT( parameters );

    auto params = std::dynamic_pointer_cast<PetscKrylovSolverParameters>( parameters );
    AMP_ASSERT( params );

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
    if ( d_bKSPCreatedInternally )
        KSPDestroy( &d_KrylovSolver );
}


/****************************************************************
 *  Initialize                                                   *
 ****************************************************************/
void PetscKrylovSolver::initialize( std::shared_ptr<SolverStrategyParameters> const params )
{
    auto parameters = std::dynamic_pointer_cast<PetscKrylovSolverParameters>( params );
    AMP_ASSERT( parameters );

    // the comm is set here of instead of the constructor because this routine
    // could be called by SNES directly
    d_comm = parameters->d_comm;
    AMP_ASSERT( !d_comm.isNull() );

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput( parameters->d_db );

    checkErr( KSPSetType( d_KrylovSolver, d_sKspType.c_str() ) );

    if ( d_KSPAppendOptionsPrefix != "" ) {
        KSPAppendOptionsPrefix( d_KrylovSolver, d_KSPAppendOptionsPrefix.c_str() );
        // PCAppendOptionsPrefix(pc, d_KSPAppendOptionsPrefix.c_str());
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
            // and the setup and apply preconditioner functions for the PCSHELL are set to
            // static member functions of this class. By doing this we do not need to introduce
            // static member functions into every SolverStrategy that might be used as a
            // preconditioner
            checkErr( PCSetType( pc, PCSHELL ) );
            checkErr( PCShellSetContext( pc, this ) );
            checkErr( PCShellSetSetUp( pc, PetscKrylovSolver::setupPreconditioner ) );
            checkErr( PCShellSetApply( pc, PetscKrylovSolver::applyPreconditioner ) );
        }
        checkErr( KSPSetPCSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
    } else {
        checkErr( PCSetType( pc, PCNONE ) );
    }

    // PetscTruth useZeroGuess = (d_bUseZeroInitialGuess) ? PETSC_TRUE : PETSC_FALSE;
    // ierr = KSPSetInitialGuessNonzero(d_KrylovSolver, useZeroGuess);

    PetscBool useNonzeroGuess = ( !d_bUseZeroInitialGuess ) ? PETSC_TRUE : PETSC_FALSE;

    checkErr( KSPSetInitialGuessNonzero( d_KrylovSolver, useNonzeroGuess ) );

    checkErr( KSPSetTolerances( d_KrylovSolver,
                                d_dRelativeTolerance,
                                d_dAbsoluteTolerance,
                                d_dDivergenceTolerance,
                                d_iMaxIterations ) );
    if ( d_bKSPCreatedInternally ) {
        checkErr( KSPSetFromOptions( d_KrylovSolver ) );
    }
    if ( d_PetscMonitor ) {
        // Add the monitor
        checkErr( KSPMonitorSet(
            d_KrylovSolver, PetscMonitor::monitorKSP, d_PetscMonitor.get(), PETSC_NULL ) );
    }
    // in this case we make the assumption we can access a PetscMat for now
    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }
}
// Function to get values from input
void PetscKrylovSolver::getFromInput( std::shared_ptr<AMP::Database> db )
{
    // fill this in
    std::string petscOptions = db->getWithDefault<std::string>( "KSPOptions", "" );
    if ( petscOptions.find( "monitor" ) != std::string::npos ) {
        petscOptions = PetscMonitor::removeMonitor( petscOptions );
        d_PetscMonitor.reset( new PetscMonitor( d_comm ) );
    }
    PetscOptionsInsertString( PETSC_NULL, petscOptions.c_str() );

    d_sKspType             = db->getWithDefault<std::string>( "ksp_type", "fgmres" );
    d_dDivergenceTolerance = db->getWithDefault<double>( "divergence_tolerance", 1.0e+03 );
    d_iMaxIterations       = db->getWithDefault<double>( "max_iterations", 1000 );

    d_KSPAppendOptionsPrefix = db->getWithDefault<std::string>( "KSPAppendOptionsPrefix", "" );

    if ( ( d_sKspType == "fgmres" ) || ( d_sKspType == "gmres" ) ) {
        d_iMaxKrylovDimension              = db->getWithDefault( "max_krylov_dimension", 20 );
        d_sGmresOrthogonalizationAlgorithm = db->getWithDefault<std::string>(
            "gmres_orthogonalization_algorithm", "modifiedgramschmidt" );
    }

    d_bUsesPreconditioner = db->getWithDefault( "uses_preconditioner", false );

    if ( d_bUsesPreconditioner ) {
        if ( db->keyExists( "pc_type" ) ) {
            d_sPcType = db->getWithDefault<std::string>( "pc_type", "none" );
        } else {
            // call error here
            AMP_ERROR( "pc_type does not exist" );
        }
        d_PcSide = db->getWithDefault<std::string>( "pc_side", "RIGHT" );
    } else {
        d_sPcType = "none";
    }
}


/****************************************************************
 *  Solve                                                        *
 ****************************************************************/
void PetscKrylovSolver::apply( std::shared_ptr<const AMP::LinearAlgebra::Vector> f,
                               std::shared_ptr<AMP::LinearAlgebra::Vector> u )
{
    PROFILE_START( "solve" );

    // Get petsc views of the vectors
    auto fVecView = AMP::LinearAlgebra::PetscVector::constView( f );
    auto uVecView = AMP::LinearAlgebra::PetscVector::view( u );

    // Check input vector states
    checkUpdateStatus( f );
    checkUpdateStatus( u );

    if ( d_iDebugPrintInfoLevel > 1 ) {
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of solution vector: " << u->L2Norm()
                  << std::endl;
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of rhs vector: " << f->L2Norm()
                  << std::endl;
    }
    Vec fVec = fVecView->getVec();
    Vec uVec = uVecView->getVec();

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
            // are set to static member functions of this class.
            // By doing this we do not need to introduce static member functions
            // into every SolverStrategy that might be used as a preconditioner
            checkErr( PCSetType( pc, PCSHELL ) );
            checkErr( PCShellSetContext( pc, this ) );
            checkErr( PCShellSetSetUp( pc, PetscKrylovSolver::setupPreconditioner ) );
            checkErr( PCShellSetApply( pc, PetscKrylovSolver::applyPreconditioner ) );
        }
        checkErr( KSPSetPCSide( d_KrylovSolver, getPCSide( d_PcSide ) ) );
    } else {
        checkErr( PCSetType( pc, PCNONE ) );
    }
    if ( d_pOperator ) {
        registerOperator( d_pOperator );
    }

    // This will replace any PETSc references to pointers we also track
    // After this, we are free to delet f_thisGetsAroundPETScSharedPtrIssue without memory leak.
    KSPSolve( d_KrylovSolver, fVec, uVec );

    if ( d_iDebugPrintInfoLevel > 2 ) {
        std::cout << "L2Norm of solution from KSP: " << u->L2Norm() << std::endl;
    }

    // Reset the solvers
    KSPReset( d_KrylovSolver );
    PROFILE_STOP( "solve" );
}


/****************************************************************
 *  Function to set the KrylovSolver                             *
 ****************************************************************/
void PetscKrylovSolver::setKrylovSolver( KSP *ksp )
{
    if ( d_bKSPCreatedInternally )
        KSPDestroy( &d_KrylovSolver );
    d_bKSPCreatedInternally = false;
    d_KrylovSolver          = *ksp;
}


/****************************************************************
 *  Function to set the register the operator                    *
 ****************************************************************/
void PetscKrylovSolver::registerOperator( const std::shared_ptr<AMP::Operator::Operator> op )
{
    // in this case we make the assumption we can access a PetscMat for now
    AMP_ASSERT( op );

    d_pOperator = op;

    auto linearOperator = std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( op );
    AMP_ASSERT( linearOperator );

    auto pMatrix = AMP::LinearAlgebra::PetscMatrix::view( linearOperator->getMatrix() );

    Mat mat;
    mat = pMatrix->getMat();

    KSPSetOperators( d_KrylovSolver, mat, mat );
}
void PetscKrylovSolver::resetOperator(
    const std::shared_ptr<AMP::Operator::OperatorParameters> params )
{
    if ( d_pOperator ) {
        d_pOperator->reset( params );
        auto linearOperator =
            std::dynamic_pointer_cast<AMP::Operator::LinearOperator>( d_pOperator );
        AMP_ASSERT( linearOperator );

        auto pMatrix = AMP::LinearAlgebra::PetscMatrix::view( linearOperator->getMatrix() );

        Mat mat;
        mat = pMatrix->getMat();

        KSPSetOperators( d_KrylovSolver, mat, mat );
    }

    // should add a mechanism for the linear operator to provide updated parameters for the
    // preconditioner operator though it's unclear where this might be necessary
    if ( d_pPreconditioner ) {
        d_pPreconditioner->resetOperator( params );
    }
}


/****************************************************************
 *  Function to setup the preconditioner                         *
 ****************************************************************/
PetscErrorCode PetscKrylovSolver::setupPreconditioner( PC pc )
{
    int ierr  = 0;
    void *ctx = nullptr;
    ierr      = PCShellGetContext( pc, &ctx );
    return ierr;
}


/****************************************************************
 *  Function to call the preconditioner                          *
 ****************************************************************/
PetscErrorCode PetscKrylovSolver::applyPreconditioner( PC pc, Vec r, Vec z )
{
    int ierr = 0;
    void *ctx;
    PCShellGetContext( pc, &ctx );
    AMP_ASSERT( ctx != nullptr );
    auto solver = reinterpret_cast<PetscKrylovSolver *>( ctx );

    auto sp_r = PETSC::getAMP( r );
    auto sp_z = PETSC::getAMP( z );

    // Make sure the vectors are in a consistent state
    sp_r->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );
    sp_z->makeConsistent( AMP::LinearAlgebra::VectorData::ScatterType::CONSISTENT_SET );

    // these tests were helpful in finding a bug
    if ( solver->getDebugPrintInfoLevel() > 5 ) {
        double norm = 0.0;
        VecNorm( r, NORM_2, &norm );
        double sp_r_norm = static_cast<double>( sp_r->L2Norm() );
        AMP_ASSERT( AMP::Utilities::approx_equal( norm, sp_r_norm ) );
    }

    // Call the preconditioner
    auto preconditioner = solver->getPreconditioner();
    if ( preconditioner != nullptr ) {
        preconditioner->apply( sp_r, sp_z );
    } else {
        // Use the identity preconditioner
        sp_z->copyVector( sp_r );
    }

    // Check for nans (no communication necessary)
    double localNorm =
        static_cast<double>( sp_z->getVectorOperations()->localL2Norm( *sp_z->getVectorData() ) );
    AMP_INSIST( localNorm == localNorm, "NaNs detected in preconditioner" );

    // not sure why, but the state of sp_z is not updated and petsc uses the cached norm
    sp_z->getVectorData()->fireDataChange();

    // these tests were helpful in finding a bug
    if ( solver->getDebugPrintInfoLevel() > 5 ) {
        double norm = 0.0;
        AMP::pout << "L2 Norm of sp_z " << sp_z->L2Norm() << std::endl;
        VecNorm( z, NORM_2, &norm );
        AMP::pout << "L2 Norm of z " << norm << std::endl;
        AMP_ASSERT( norm == sp_z->L2Norm() );
    }

    return ( ierr );
}


} // namespace Solver
} // namespace AMP
