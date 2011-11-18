#include <iostream>
#include "LinearTimeOperator.h"

#include "utils/Utilities.h"

#ifndef included_TimeIntegratorParameters
#include "TimeIntegratorParameters.h"
#endif

#ifndef included_IDATimeIntegrator
#include "IDATimeIntegrator.h"
#endif

#ifndef included_AMP_OperatorBuilder
#include "operators/OperatorBuilder.h"
#endif

#ifdef USE_SUNDIALS
extern "C"{
#include "ida/ida_impl.h"
}

namespace AMP{
namespace TimeIntegrator{
    
    
#define AMPVEC_CAST(v) (static_cast<ManagedSundialsVector*>(v->content))    
    
    /*
     ************************************************************************
     *                                                                      *
     *  Constructor.                                                        *
     *                                                                      *
     ************************************************************************
     */
    IDATimeIntegrator::IDATimeIntegrator( boost::shared_ptr< TimeIntegratorParameters > parameters ):TimeIntegrator(parameters)
    {
      initialize( parameters );
    }
    
    /*
     ************************************************************************
     *                                                                      *
     *  Destructor.                                                         *
     *                                                                      *
     ************************************************************************
     */
    IDATimeIntegrator::~IDATimeIntegrator()
    {
    }
    
    /*
     ************************************************************************
     *                                                                      *
     * Initialize.                                                          *
     *                                                                      *
     ************************************************************************
     */
    void
    IDATimeIntegrator::initialize( boost::shared_ptr< TimeIntegratorParameters> parameters )
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(parameters.get() != NULL);
#endif
        getFromInput( parameters->d_db );
        
        boost::shared_ptr< IDATimeIntegratorParameters > params = boost::dynamic_pointer_cast< IDATimeIntegratorParameters > (parameters);
        d_solution_prime = (params->d_ic_vector_prime)->cloneVector();
        d_solution_prime->copyVector(params->d_ic_vector_prime);
        
        d_pPreconditioner = params->d_pPreconditioner;

        // reuse the time integrator database, and put additional fields in
        boost::shared_ptr<AMP::Database> timeOperator_db = params->d_db;
        timeOperator_db->putDouble("CurrentDt", d_current_dt);
        timeOperator_db->putDouble("CurrentTime", d_current_time);
        timeOperator_db->putString("name", "TimeOperator");

        // setup the parameter object for the IDATimeOperator
        boost::shared_ptr<AMP::TimeIntegrator::IDATimeOperatorParameters> idaTimeOp_Params ( new AMP::TimeIntegrator::IDATimeOperatorParameters(params->d_db));
        
        idaTimeOp_Params->d_pRhsOperator = parameters->d_operator;
        idaTimeOp_Params->d_pMassOperator = parameters->d_pMassOperator;
        idaTimeOp_Params->d_pSourceTerm = parameters->d_pSourceTerm;
        idaTimeOp_Params->d_pAlgebraicVariable = parameters->d_pAlgebraicVariable;
          
        // create the time operator
        boost::shared_ptr<AMP::TimeIntegrator::IDATimeOperator> idaTimeOp(new AMP::TimeIntegrator::IDATimeOperator(idaTimeOp_Params));
        d_pIDATimeOperator = idaTimeOp;
        
        // if we want to create the LinearTimeOperator internally
        if(d_createLinearOperatorInternally)
          {
            if(d_bLinearRhsOperator && d_bLinearMassOperator)
              {
            boost::shared_ptr<AMP::TimeIntegrator::TimeOperatorParameters> linearTimeOperatorParams = boost::dynamic_pointer_cast<AMP::TimeIntegrator::TimeOperatorParameters>(idaTimeOp->getJacobianParameters(d_solution));
                        boost::shared_ptr<AMP::Database> timeOperator_db = linearTimeOperatorParams->d_db;
                        timeOperator_db->putDouble("CurrentDt", d_current_dt);
                        timeOperator_db->putDouble("CurrentTime", d_current_time);
                        timeOperator_db->putString("name", "TimeOperator");

            d_pLinearTimeOperator.reset(new LinearTimeOperator(linearTimeOperatorParams));
              }
            else
              {
            AMP::pout << "ERROR: IDATimeIntegrator::initialize(): creation of linear time operators internally is only currently supported for linear mass and rhs operators" << std::endl;
              }
            
            AMP_INSIST(d_pPreconditioner.get()!=NULL, "ERROR: IDATimeIntegrator::initialize(): creation of linear time operators internally is only currently supported with a valid non NULL preconditioner ");
    
            d_pPreconditioner->registerOperator(d_pLinearTimeOperator);
            AMP::pout << " linear op being created internally" << std::endl;            
          }
        /*
         * Initialize data members from input.
         */
        
        setupVectors();    
        initializeIDA();
        
    }
    
    void IDATimeIntegrator::initializeIDA()
    {
        N_Vector id=NULL;
        
        AMP::LinearAlgebra::Vector::shared_ptr  pSundials_sol = AMP::LinearAlgebra::SundialsVector::view ( d_solution );
        AMP::LinearAlgebra::Vector::shared_ptr  pSundials_sol_prime = AMP::LinearAlgebra::SundialsVector::view ( d_solution_prime );
        
        
        id = N_VClone(pSundials_sol->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector());
        
        d_ida_mem = IDACreate();
        assert(d_ida_mem!=0);    
        
        int ierr = IDASetUserData(d_ida_mem, this);
        assert(ierr==IDA_SUCCESS);
        
        N_VConst(1.0, id);
        ierr = IDASetId(d_ida_mem, id);
        assert(ierr==IDA_SUCCESS);
        
        // boost::shared_ptr<AMP::LinearAlgebra::SundialsVector> pSun_nvec =  boost::dynamic_pointer_cast<AMP::LinearAlgebra::SundialsVector>(pSundials_sol);
                
        ierr = IDAInit(d_ida_mem, IDAResTrial, d_initial_time, pSundials_sol->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector(), pSundials_sol_prime->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector());
        assert(ierr==IDA_SUCCESS);        
        
        ierr = IDASStolerances(d_ida_mem, d_relative_tolerance, d_absolute_tolerance);
        assert(ierr==IDA_SUCCESS);
                
        N_VDestroy(id);
        
        // set the initial step size
        ierr = IDASetInitStep(d_ida_mem, d_initial_dt);
        assert(ierr==IDA_SUCCESS);

        // set the max step size
        ierr = IDASetMaxStep(d_ida_mem, d_max_dt);
        assert(ierr==IDA_SUCCESS);

        // set final time
        ierr = IDASetStopTime(d_ida_mem, d_final_time);
        assert(ierr==IDA_SUCCESS);
        
        // ideally, the linear solver type needs to be determined by the user input
        ierr = IDASpgmr(d_ida_mem, 0);
        assert(ierr==IDA_SUCCESS);

        if(d_bUsePreconditioner)
          {
            ierr = IDASpilsSetPreconditioner(d_ida_mem, IDAPrecSetup, IDAPrecSolve);
            assert(ierr==IDASPILS_SUCCESS);
          }
        
        //ierr = IDASpilsSetMaxRestarts(d_ida_mem, 15);
        ierr = IDASpilsSetMaxRestarts(d_ida_mem, 100);
        assert(ierr==IDA_SUCCESS);


                

        // if we want IDA to calculate consistent IC's
        
        if(d_bCallCalcIC)

          {
            double tout1 = d_initial_time+d_current_dt;
            // choice of the last argument...?
            ierr = IDACalcIC(d_ida_mem, IDA_YA_YDP_INIT, tout1);    
            assert(ierr==IDA_SUCCESS);
            ierr = IDAGetConsistentIC(d_ida_mem, pSundials_sol->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector(), pSundials_sol_prime->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector());
            assert(ierr==IDA_SUCCESS);
        
          }
        
        // For now, just use BDF2
        //ierr = IDASetMaxOrd(d_ida_mem, 2);
        //assert(ierr==IDA_SUCCESS);

        
        //d_init_step_size = 0.1;
        //ierr = IDASetInitStep(d_ida_mem, d_init_step_size);
        //assert(ierr==IDA_SUCCESS);
        
    }
    
    void
    IDATimeIntegrator::reset( boost::shared_ptr< TimeIntegratorParameters > parameters )
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(parameters.get() != NULL);
#endif
        
        abort();
    }
    
    void IDATimeIntegrator::setupVectors( void )
    {
        
        // clone vectors so they have the same data layout as d_solution 
        d_residual   = d_solution->cloneVector();
        /*
         * Set initial value of vectors to 0.
         */
        d_residual->setToScalar((double) 0.0);
        
    }
    
    /*
     ************************************************************************
     *                                                                      *
     *  Update internal state to reflect time advanced solution.            *
     *                                                                      *
     ************************************************************************
     */
    void
    IDATimeIntegrator::updateSolution( void )
    {
        /*
         int retval = IDA_SUCCESS;
         realtype hlast;
         
         retval = IDAGetLastStep(d_ida_mem, &hlast);
         assert(retval==IDA_SUCCESS);
         
         d_current_time += hlast;
         //d_solution->add( d_predictor, d_corrector );
         */
    }
    
    
    /*
     ************************************************************************
     *                                                                      *
     * Read input from database.                                            *
     *                                                                      *
     ************************************************************************
     */
    void
    IDATimeIntegrator::getFromInput( boost::shared_ptr<AMP::Database> input_db )
    {
        if ( input_db->keyExists("bLinearMassOperator") ) {
            d_bLinearMassOperator = input_db->getBool("bLinearMassOperator");
        } else {
            AMP_ERROR(d_object_name << " -- Key data `bLinearMassOperator'"
                       << " missing in input.");
        }
        
        if ( input_db->keyExists("bLinearRhsOperator") ) {
            d_bLinearRhsOperator = input_db->getBool("bLinearRhsOperator");
        } else {
            AMP_ERROR(d_object_name << " -- Key data `bLinearRhsOperator'"
                       << " missing in input.");
        }
        
        if ( input_db->keyExists("linear_solver_type") ) {
            d_linear_solver_type = input_db->getInteger("linear_solver_type");
        } else {
            AMP_ERROR(d_object_name << " -- Key data `linear_solver_type'"
                       << " missing in input.");
        }
        
        if ( input_db->keyExists("relative_tolerance") ) {
            d_relative_tolerance = input_db->getDouble("relative_tolerance");
        } else {
            AMP_ERROR(d_object_name << " -- Key data `relative_tolerance'"
                       << " missing in input.");
        }
        
        if ( input_db->keyExists("absolute_tolerance") ) {
            d_absolute_tolerance = input_db->getDouble("absolute_tolerance");
        } else {
            AMP_ERROR(d_object_name << " -- Key data `absolute_tolerance'"
                       << " missing in input.");
        }

        d_bCallCalcIC = input_db->getBoolWithDefault("CallCalcIC", true);
        d_bUsePreconditioner = input_db->getBoolWithDefault("usePreconditioner", true);
        
        d_createLinearOperatorInternally = input_db->getBoolWithDefault("createLinearTimeOperatorInternally", false);
        d_bManufacturedProblem = input_db->getBoolWithDefault("bManufacturedProblem", false);
    }
    
    double 
    IDATimeIntegrator::getNextDt(const bool good_solution)
    {
        int ierr = IDAGetCurrentStep(d_ida_mem, &d_current_dt);
        assert(ierr!=IDA_SUCCESS);
        
        return d_current_dt;
    }
    
    int IDATimeIntegrator::advanceSolution( const double dt, const bool first_step )
    {
        int retval = IDA_SUCCESS;
        double hin_actual=0.0;
        double hlast=0.0;
        double hcur=0.0;
        long int nniters=0;
        long int nncfails=0;
        long int netfails=0;
        long int nliters=0;
        long int npsolves=0;
        
        AMP::LinearAlgebra::Vector::shared_ptr  ptr_y = AMP::LinearAlgebra::SundialsVector::view( d_solution );
        AMP::LinearAlgebra::Vector::shared_ptr  ptr_ydot = AMP::LinearAlgebra::SundialsVector::view( d_solution_prime );

        // specify some initial step
        hcur=d_current_time+dt;

        AMP::pout << "before IDASolve" << std::endl;
        retval = IDASolve(d_ida_mem, hcur, &d_current_time, ptr_y->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector(), ptr_ydot->castTo<AMP::LinearAlgebra::SundialsVector>().getNVector(), IDA_ONE_STEP);
        //should be fixed.
        AMP::pout << "after IDASolve" << std::endl;
        retval = IDAGetActualInitStep(d_ida_mem, &hin_actual);
        //cout << "hin_actual = " << hin_actual << endl;
        
        retval = IDAGetLastStep(d_ida_mem, &hlast);
        AMP::pout << "hlast = " << hlast << std::endl;
        
        retval = IDAGetCurrentStep(d_ida_mem, &hcur);
        d_current_dt = hcur;
        AMP::pout << "hcur = " << hcur << std::endl;
        
        retval = IDAGetNumNonlinSolvIters(d_ida_mem, &nniters);
        AMP::pout << "nniters = " << nniters << std::endl;
        
        retval = IDAGetNumNonlinSolvConvFails(d_ida_mem, &nncfails);
        AMP::pout << "nncfails = " << nncfails << std::endl;
        
        retval = IDAGetNumErrTestFails(d_ida_mem, &netfails);
        AMP::pout << "netfails = " << netfails << std::endl;
        
        retval = IDASpilsGetNumPrecSolves(d_ida_mem, &npsolves);
        AMP::pout << "npsolves = " << npsolves << std::endl;
        
        retval = IDASpilsGetNumLinIters(d_ida_mem, &nliters);
        AMP::pout << "nliters = " << nliters << std::endl;
        //cout << "current time = " << d_current_time << endl;
        
        //double maxval;
        //cout << "maxval = " << ptr_y->getNVector()->ops->nvmaxnorm(ptr_y->getNVector()) << endl;
        
        // solution has been just updated; the source term should be updated accordingly.    
        //updateSourceTerm();
        
        return(retval);
    }
    
    void IDATimeIntegrator::updateSourceTerm( void ) 
    {
        // Ideally, IDATimeIntegrator should call a function which takes d_solution and t and returns the appropriate 
        // updated source term.
        
    }
    /*
     ************************************************************************
     *                                                                      *
     *  Check whether time advanced solution is acceptable.                 *
     *                                                                      *
     ************************************************************************
     */
    bool
    IDATimeIntegrator::checkNewSolution( void ) const
    {
        /*
         * Ordinarily we would check the actual error in the solution
         * (proportional to the size of d_corrector) against a specified
         * tolerance.  For now, accept everything.
         */
        return(true);
    }
    
    
    int IDATimeIntegrator::IDAResTrial(realtype tt, N_Vector yy, N_Vector yp,
                                       N_Vector rr, void *user_data)
    {
        
        boost::shared_ptr<AMP::LinearAlgebra::Vector> f;
        
        
        AMP::LinearAlgebra::ExternalVectorDeleter d;


        AMP::LinearAlgebra::Vector * pyy = static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(yy->content);
        AMP::LinearAlgebra::Vector * pyp = static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(yp->content);
        boost::shared_ptr<AMP::LinearAlgebra::Vector> amp_yy (pyy, d);
        boost::shared_ptr<AMP::LinearAlgebra::Vector> amp_yp (pyp, d);
        
        ((IDATimeIntegrator*)user_data)->getIDATimeOperator()->registerIDATimeDerivative(amp_yp);
        
        AMP::LinearAlgebra::Vector::shared_ptr d_residual_Sundials = ((IDATimeIntegrator *)user_data)->getResidualVector();

                double currentTime = ((IDATimeIntegrator *)user_data)->getCurrentTime();
        ((IDATimeIntegrator*)user_data)->getIDATimeOperator()->registerCurrentTime(currentTime);
    
        ((IDATimeIntegrator*)user_data)->getIDATimeOperator()->apply(f, amp_yy, d_residual_Sundials, 1.0, 0.0);

        (static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(rr->content))->copyVector( *d_residual_Sundials );
        
        
        return(IDA_SUCCESS);
        
    }
    
    int IDATimeIntegrator::IDAPrecSetup(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
                                        realtype cj, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {
        int last_order, current_order;
        double last_stepsize, current_stepsize;
        
        int ierr = IDAGetLastOrder(((IDATimeIntegrator*)user_data)->getIDAMem(), &last_order);
        ierr = IDAGetCurrentOrder(((IDATimeIntegrator*)user_data)->getIDAMem(), &current_order);
        ierr = IDAGetLastStep(((IDATimeIntegrator*)user_data)->getIDAMem(), &last_stepsize);
        ierr = IDAGetCurrentStep(((IDATimeIntegrator*)user_data)->getIDAMem(), &current_stepsize);
        assert(ierr==IDA_SUCCESS);
        
        AMP::pout << "cj = " << cj << std::endl;
        
        
        
        bool OrderHasChanged = FALSE, StepSizeHasChanged = FALSE;
        if (last_order != current_order){
            OrderHasChanged = TRUE;
        }
        if (last_stepsize != current_stepsize){
            StepSizeHasChanged = TRUE;
        }
        
        
        if (OrderHasChanged || StepSizeHasChanged)
        {
          
          AMP::LinearAlgebra::ExternalVectorDeleter d;

          AMP::LinearAlgebra::Vector * pyy = static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(yy->content);
          boost::shared_ptr<AMP::LinearAlgebra::Vector> amp_yy (pyy, d);
          boost::shared_ptr<AMP::Operator::OperatorParameters> jacParams = ((IDATimeIntegrator*)user_data)->getIDATimeOperator()->getJacobianParameters(amp_yy);
          boost::shared_ptr<AMP::Database> &db = jacParams->d_db;
          db->putDouble("ScalingFactor", cj);
          db->putDouble("CurrentTime", tt);
          boost::shared_ptr<AMP::Solver::SolverStrategy> pSolver = ((IDATimeIntegrator*)user_data)->getPreconditioner();
          //double currentTime = ((IDATimeIntegrator *)user_data)->getCurrentTime();

          // BP:  commenting out these lines that appear to be buggy
          // I will need to fix these lines so that the current time is set through the database object
          // we can have column time operators and in that case the lines below fail
          //          if (((IDATimeIntegrator *)user_data)->getBoolManufacturedProblem()) {
          //          boost::shared_ptr<AMP::TimeIntegrator::LinearTimeOperator> pLinearTimeOperator = boost::dynamic_pointer_cast<AMP::TimeIntegrator::LinearTimeOperator> ( pSolver->getOperator() );          
          //          AMP_INSIST(pLinearTimeOperator!=NULL, "In IDATimeIntegrator::IDAPrecSetup - pLinearTimeOperator is NULL");
          //          pLinearTimeOperator->registerCurrentTime(currentTime);
          pSolver->resetOperator(jacParams);
          //          }
          
        }
        
        return(IDA_SUCCESS);
    }
    
    
    int IDATimeIntegrator::IDAPrecSolve(realtype tt, N_Vector yy, N_Vector yp, N_Vector rr,
                                        N_Vector rvec, N_Vector zvec,
                                        realtype cj, realtype delta,
                                        void *user_data, N_Vector tmp)
    {
        AMP::LinearAlgebra::ExternalVectorDeleter d;
        AMP::LinearAlgebra::Vector * pr = static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(rvec->content);
        AMP::LinearAlgebra::Vector * pz = static_cast<AMP::LinearAlgebra::ManagedSundialsVector*>(zvec->content);
        boost::shared_ptr<AMP::LinearAlgebra::Vector> amp_rvec (pr, d);
        boost::shared_ptr<AMP::LinearAlgebra::Vector> amp_zvec (pz, d);
        
        ((IDATimeIntegrator*)user_data)->getPreconditioner()->solve(amp_rvec, amp_zvec);
        
        return(IDA_SUCCESS);
    }
    
    boost::shared_ptr<IDATimeOperator> IDATimeIntegrator::getIDATimeOperator() const
    {
        return(d_pIDATimeOperator);
    }
    
    boost::shared_ptr<LinearTimeOperator> IDATimeIntegrator::getLinearTimeOperator() const
    {
        return(d_pLinearTimeOperator);
    }
    
    
    boost::shared_ptr<AMP::LinearAlgebra::Vector> IDATimeIntegrator::getResidualVector() const
    {
        return(d_residual);
        
    }
    
    
    boost::shared_ptr<AMP::LinearAlgebra::Vector> IDATimeIntegrator::getSourceTerm() const
    {
        return(d_pSourceTerm);
    }
    
    
    
    
}    
}

#endif







