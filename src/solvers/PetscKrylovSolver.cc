#include "PetscKrylovSolver.h"
#include "vectors/petsc/ManagedPetscVector.h"
#include "operators/LinearOperator.h"
#include "matrices/Matrix.h"
#include "matrices/petsc/PetscMatrix.h"
#include "utils/ProfilerApp.h"

extern "C"{
#include "assert.h"
#include "petscpc.h"
}

namespace AMP {
namespace Solver {


#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
static inline void checkErr(int ierr) {
    AMP_INSIST(ierr==0, "Petsc returned non-zero error code");
}
#elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
static inline void checkErr(PetscErrorCode ierr) {
    AMP_INSIST(ierr==0, "Petsc returned non-zero error code");
}
#else
    #error Not programmed for this version yet
#endif



/****************************************************************
*  Constructors                                                 *
****************************************************************/
PetscKrylovSolver::PetscKrylovSolver()
{
    d_bKSPCreatedInternally = false;
    d_KrylovSolver = NULL;
}
PetscKrylovSolver::PetscKrylovSolver(boost::shared_ptr<PetscKrylovSolverParameters> parameters):SolverStrategy(parameters)
{
    assert(parameters.get()!=NULL);

    // Create a default KrylovSolver
    d_bKSPCreatedInternally = true;
    KSPCreate( parameters->d_comm.getCommunicator(), &d_KrylovSolver );

    // Initialize
    initialize(parameters);
}


/****************************************************************
*  De-constructor                                               *
****************************************************************/
PetscKrylovSolver::~PetscKrylovSolver()
{
    if(d_bKSPCreatedInternally) {
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            KSPDestroy(d_KrylovSolver);
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            KSPDestroy(&d_KrylovSolver);
        #else
            #error Not programmed for this version yet
        #endif
    }
}


/****************************************************************
*  Initialize                                                   *
****************************************************************/
void
PetscKrylovSolver::initialize(boost::shared_ptr<SolverStrategyParameters> const params)
{
    boost::shared_ptr<PetscKrylovSolverParameters> parameters = boost::dynamic_pointer_cast<PetscKrylovSolverParameters>(params);
    AMP_ASSERT(parameters.get()!=NULL);
    d_comm = parameters->d_comm;
    AMP_ASSERT(!d_comm.isNull());

    d_pPreconditioner = parameters->d_pPreconditioner;

    getFromInput(parameters->d_db);

    checkErr(KSPSetType(d_KrylovSolver,d_sKspType.c_str()));

    if(d_KSPAppendOptionsPrefix!="") {
        KSPAppendOptionsPrefix(d_KrylovSolver, d_KSPAppendOptionsPrefix.c_str());
        //      PCAppendOptionsPrefix(pc, d_KSPAppendOptionsPrefix.c_str());
    }

    if((d_sKspType=="fgmres")||(d_sKspType=="gmres")) {
        checkErr(KSPGMRESSetRestart(d_KrylovSolver, d_iMaxKrylovDimension));
    }

    // Create the preconditioner
    PC pc;
    checkErr(KSPGetPC(d_KrylovSolver,&pc));
    if(d_bUsesPreconditioner) {
        if(d_sPcType !="shell") {
            // the pointer to the preconditioner should be NULL if we are using a Petsc internal PC
            assert(d_pPreconditioner.get()==NULL);
            PCSetType(pc, d_sPcType.c_str());
        } else {
            // for a shell preconditioner the user context is set to an instance of this class
            // and the setup and apply preconditioner functions for the PCSHELL
            // are set to static member functions of this class. By doing this we do not need to introduce
            // static member functions into every SolverStrategy that might be used as a preconditioner
            checkErr(PCSetType(pc,PCSHELL));
            checkErr(PCShellSetContext(pc, this));
            checkErr(PCShellSetSetUp(pc, PetscKrylovSolver::setupPreconditioner));
            checkErr(PCShellSetApply(pc, PetscKrylovSolver::applyPreconditioner));
        }
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            checkErr(KSPSetPreconditionerSide(d_KrylovSolver, d_PcSide));
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            checkErr(KSPSetPCSide(d_KrylovSolver, d_PcSide));
        #else
            #error Not programmed for this version yet
        #endif
    } else {
        checkErr(PCSetType(pc,PCNONE));
    }

    //PetscTruth useZeroGuess = (d_bUseZeroInitialGuess) ? PETSC_TRUE : PETSC_FALSE;
    //ierr = KSPSetInitialGuessNonzero(d_KrylovSolver, useZeroGuess);

    #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
        PetscTruth useNonzeroGuess = (!d_bUseZeroInitialGuess) ? PETSC_TRUE : PETSC_FALSE;
    #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
        PetscBool useNonzeroGuess = (!d_bUseZeroInitialGuess) ? PETSC_TRUE : PETSC_FALSE;
    #else
        #error Not programmed for this version yet
    #endif
    checkErr(KSPSetInitialGuessNonzero(d_KrylovSolver, useNonzeroGuess));

    checkErr(KSPSetTolerances(d_KrylovSolver, d_dRelativeTolerance, d_dAbsoluteTolerance, d_dDivergenceTolerance, d_iMaxIterations));
    if(d_bKSPCreatedInternally){
        checkErr(KSPSetFromOptions(d_KrylovSolver));
    }
    if ( d_PetscMonitor.get()!=NULL ) {
        // Add the monitor
        checkErr(KSPMonitorSet(d_KrylovSolver,PetscMonitor::monitorKSP,d_PetscMonitor.get(),PETSC_NULL));
    }
    // in this case we make the assumption we can access a PetscMat for now
    if(d_pOperator.get()!=NULL) {
       registerOperator(d_pOperator);
    }
}
// Function to get values from input
void PetscKrylovSolver::getFromInput(const boost::shared_ptr<AMP::Database> &db)
{
    // fill this in
    std::string petscOptions = db->getStringWithDefault("KSPOptions", "");
    if ( petscOptions.find("monitor")!=std::string::npos ) {
        petscOptions = PetscMonitor::removeMonitor(petscOptions);
        d_PetscMonitor.reset( new PetscMonitor( d_comm ) );
    }
    PetscOptionsInsertString(petscOptions.c_str());

    d_sKspType = db->getStringWithDefault("ksp_type", "fgmres");
    d_dRelativeTolerance = db->getDoubleWithDefault("relative_tolerance", 1.0e-9);
    d_dAbsoluteTolerance = db->getDoubleWithDefault("absolute_tolerance", 1.0e-14);
    d_dDivergenceTolerance = db->getDoubleWithDefault("divergence_tolerance", 1.0e+03);
    d_iMaxIterations       = db->getDoubleWithDefault("max_iterations", 1000);

    d_KSPAppendOptionsPrefix = db->getStringWithDefault("KSPAppendOptionsPrefix", "");

    if((d_sKspType=="fgmres")||(d_sKspType=="gmres")) {
        d_iMaxKrylovDimension = db->getIntegerWithDefault("max_krylov_dimension", 20);
        d_sGmresOrthogonalizationAlgorithm = db->getStringWithDefault("gmres_orthogonalization_algorithm", "modifiedgramschmidt");
    }

    d_bUsesPreconditioner = db->getBoolWithDefault("uses_preconditioner", false);

    if (d_bUsesPreconditioner) {
        if (db->keyExists("pc_type")) {
            d_sPcType = db->getStringWithDefault("pc_type", "none");
        } else {
            // call error here
            AMP_ERROR("pc_type does not exist");
        }
        std::string pc_side = db->getStringWithDefault("pc_side", "RIGHT");
        if (pc_side=="RIGHT") {
            d_PcSide = PC_RIGHT;
        } else if (pc_side=="LEFT") {
            d_PcSide = PC_LEFT;
        } else if (pc_side=="SYMMETRIC") {
            d_PcSide = PC_SYMMETRIC;
        } else {
            AMP_ERROR("Unknown value for pc_type");
        }
    } else {
        d_sPcType = "none";
    }

}

/****************************************************************
*  Solve                                                        *
****************************************************************/
void
PetscKrylovSolver::solve(boost::shared_ptr<const AMP::LinearAlgebra::Vector>  f,
                  boost::shared_ptr<AMP::LinearAlgebra::Vector>  u)
{
    PROFILE_START("solve");
    #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
        // fVecView and uVecView may be held in KSPSolve internals.
        // by declaring a temporary vector, we ensure that the KSPSolve
        // will be replaced by fVecView and uVecView before they are
        // destroyed by boost.
        AMP::LinearAlgebra::Vector::const_shared_ptr  f_thisGetsAroundPETScSharedPtrIssue = fVecView;
        AMP::LinearAlgebra::Vector::shared_ptr  u_thisGetsAroundPETScSharedPtrIssue = uVecView;
    #endif

    // Get petsc views of the vectors
    #if !( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
        AMP::LinearAlgebra::Vector::const_shared_ptr  fVecView;
        AMP::LinearAlgebra::Vector::shared_ptr  uVecView;
    #endif
    fVecView = AMP::LinearAlgebra::PetscVector::constView ( f );
    uVecView = AMP::LinearAlgebra::PetscVector::view ( u );

    // Check input vector states
    AMP_ASSERT( (f->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (f->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    AMP_ASSERT( (u->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (u->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    AMP_ASSERT( (fVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (fVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    AMP_ASSERT( (uVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (uVecView->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );

    if(d_iDebugPrintInfoLevel>1) {
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of solution vector: " << u->L2Norm() << std::endl;
        std::cout << "PetscKrylovSolver::solve: initial L2Norm of rhs vector: " << f->L2Norm() << std::endl;
    }
    Vec fVec = fVecView->castTo<AMP::LinearAlgebra::PetscVector>().getVec();
    Vec uVec = uVecView->castTo<AMP::LinearAlgebra::PetscVector>().getVec();

    // Create the preconditioner and re-register the operator
    PC pc;
    checkErr(KSPGetPC(d_KrylovSolver,&pc));
    if(d_bUsesPreconditioner) {
        if(d_sPcType !="shell") {
            // the pointer to the preconditioner should be NULL if we are using a Petsc internal PC
            assert(d_pPreconditioner.get()==NULL);
            PCSetType(pc, d_sPcType.c_str());
        } else {
            // for a shell preconditioner the user context is set to an instance of this class
            // and the setup and apply preconditioner functions for the PCSHELL
            // are set to static member functions of this class. By doing this we do not need to introduce
            // static member functions into every SolverStrategy that might be used as a preconditioner
            checkErr(PCSetType(pc,PCSHELL));
            checkErr(PCShellSetContext(pc, this));
            checkErr(PCShellSetSetUp(pc, PetscKrylovSolver::setupPreconditioner));
            checkErr(PCShellSetApply(pc, PetscKrylovSolver::applyPreconditioner));
        }
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            checkErr(KSPSetPreconditionerSide(d_KrylovSolver, d_PcSide));
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            checkErr(KSPSetPCSide(d_KrylovSolver, d_PcSide));
        #else
            #error Not programmed for this version yet
        #endif
    } else {
        checkErr(PCSetType(pc,PCNONE));
    }
    if(d_pOperator.get()!=NULL) {
       registerOperator(d_pOperator);
    }

    // This will replace any PETSc references to pointers we also track
    // After this, we are free to delet f_thisGetsAroundPETScSharedPtrIssue without memory leak.
    KSPSolve(d_KrylovSolver, fVec, uVec);

    if(d_iDebugPrintInfoLevel>2) {
        std::cout << "L2Norm of solution from KSP: " << u->L2Norm() << std::endl;
    }

    // Reset the solvers
    #if !( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
        KSPReset(d_KrylovSolver);
    #endif
    PROFILE_STOP("solve");
}


/****************************************************************
*  Function to set the KrylovSolver                             *
****************************************************************/
void PetscKrylovSolver::setKrylovSolver(KSP *ksp)
{
    if (d_bKSPCreatedInternally) {
        #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
            KSPDestroy(d_KrylovSolver);
        #elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
            KSPDestroy(&d_KrylovSolver);
        #else
            #error Not programmed for this version yet
        #endif
    }
    d_bKSPCreatedInternally = false;
    d_KrylovSolver = *ksp;
}


/****************************************************************
*  Function to set the register the operator                    *
****************************************************************/
void PetscKrylovSolver::registerOperator(const boost::shared_ptr<AMP::Operator::Operator> op)
{
  // in this case we make the assumption we can access a PetscMat for now
  assert(op.get()!=NULL);

  d_pOperator = op;

  boost::shared_ptr<AMP::Operator::LinearOperator> linearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(op);
  assert(linearOperator.get() != NULL);

  boost::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix = boost::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(linearOperator->getMatrix());
  assert(pMatrix.get()!=NULL);

  Mat mat;
  mat = pMatrix->getMat();

  KSPSetOperators(d_KrylovSolver,mat,mat,DIFFERENT_NONZERO_PATTERN);

}
void PetscKrylovSolver::resetOperator(const boost::shared_ptr<AMP::Operator::OperatorParameters> params)
{
  if(d_pOperator.get()!=NULL) {
      d_pOperator->reset(params);
      boost::shared_ptr<AMP::Operator::LinearOperator> linearOperator = boost::dynamic_pointer_cast<AMP::Operator::LinearOperator>(d_pOperator);
      assert(linearOperator.get() != NULL);

      boost::shared_ptr<AMP::LinearAlgebra::PetscMatrix> pMatrix = boost::dynamic_pointer_cast<AMP::LinearAlgebra::PetscMatrix>(linearOperator->getMatrix());
      assert(pMatrix.get()!=NULL);

      Mat mat;
      mat = pMatrix->getMat();

      // for now we will assume the pc is affected at every iteration
      KSPSetOperators(d_KrylovSolver,mat,mat,DIFFERENT_NONZERO_PATTERN);
  }

  // should add a mechanism for the linear operator to provide updated parameters for the preconditioner operator
  // though it's unclear where this might be necessary
  if(d_pPreconditioner.get()!=NULL) {
      d_pPreconditioner->resetOperator(params);
  }
}


/****************************************************************
*  Function to setup the preconditioner                         *
****************************************************************/
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
int PetscKrylovSolver::setupPreconditioner(void*)
{
   int ierr = 0;
   //return( ((PetscKrylovSolver*)ctx)->getPreconditioner()->reset() );
   return ierr;
}
#elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
PetscErrorCode PetscKrylovSolver::setupPreconditioner(PC pc)
{
   int ierr = 0;
   void *ctx = NULL;
   ierr = PCShellGetContext(pc, &ctx);
   return ierr;
}
#else
    #error Not programmed for this version yet
#endif


/****************************************************************
*  Function to call the preconditioner                          *
****************************************************************/
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==0 )
PetscErrorCode  PetscKrylovSolver::applyPreconditioner(void* ctx, Vec r, Vec z)
#elif ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
PetscErrorCode  PetscKrylovSolver::applyPreconditioner(PC pc, Vec r, Vec z)
#else
    #error Not programmed for this version yet
#endif
{
    int ierr = 0;
    #if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==2 )
        void* ctx;
        PCShellGetContext(pc,&ctx);
    #endif
    AMP_ASSERT(ctx!=NULL);

    boost::shared_ptr<AMP::LinearAlgebra::Vector> sp_r ( reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>(r->data) , AMP::LinearAlgebra::ExternalVectorDeleter() );
    boost::shared_ptr<AMP::LinearAlgebra::Vector> sp_z ( reinterpret_cast<AMP::LinearAlgebra::ManagedPetscVector *>(z->data) , AMP::LinearAlgebra::ExternalVectorDeleter() );

    // these tests were helpful in finding a bug
    if(((PetscKrylovSolver*)ctx)->getDebugPrintInfoLevel()>5) {
        double norm=0.0;
        VecNorm(r, NORM_2,&norm);
        double sp_r_norm = sp_r->L2Norm();
        AMP_ASSERT(AMP::Utilities::approx_equal(norm, sp_r_norm));
    }  
  
    AMP_ASSERT( (sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (sp_r->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );
    AMP_ASSERT( (sp_z->getUpdateStatus() == AMP::LinearAlgebra::Vector::UNCHANGED) ||
        (sp_z->getUpdateStatus() == AMP::LinearAlgebra::Vector::LOCAL_CHANGED) );

    // Call the preconditioner
    boost::shared_ptr<AMP::Solver::SolverStrategy> preconditioner = 
        ((PetscKrylovSolver*)ctx)->getPreconditioner();
    if ( preconditioner!=NULL ) {
        preconditioner->solve(sp_r,sp_z);
    } else {
        // Use the identity preconditioner
        sp_z->copyVector(sp_r);
    }

    // Check for nans (no communication necessary)
    double localNorm = sp_z->localL2Norm();
    AMP_INSIST(localNorm==localNorm,"NaNs detected in preconditioner");

    // not sure why, but the state of sp_z is not updated
    // and petsc uses the cached norm
    if ( sp_z->isA<AMP::LinearAlgebra::DataChangeFirer>() ) {
        sp_z->castTo<AMP::LinearAlgebra::DataChangeFirer>().fireDataChange();
    }

    // these tests were helpful in finding a bug
    if(((PetscKrylovSolver*)ctx)->getDebugPrintInfoLevel()>5) {
        double norm=0.0;
        AMP::pout << "L2 Norm of sp_z " << sp_z->L2Norm() << std::endl;
        VecNorm(z, NORM_2,&norm);
        AMP::pout << "L2 Norm of z " << norm << std::endl;
        AMP_ASSERT(norm==sp_z->L2Norm());      
    }
  
    return (ierr);
}



}
}

