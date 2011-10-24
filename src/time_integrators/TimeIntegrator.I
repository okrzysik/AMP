#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace AMP{
namespace TimeIntegrator{

inline
double TimeIntegrator::getCurrentDt() const
{
   return(d_current_dt);
}

inline
double TimeIntegrator::getInitialTime() const
{
   return(d_initial_time);
}

inline
double TimeIntegrator::getFinalTime() const
{
   return(d_final_time);
}

inline
double TimeIntegrator::getCurrentTime() const
{
   return (d_current_time);
}

inline
bool TimeIntegrator::stepsRemaining() const
{
   return(d_integrator_step < d_max_integrator_steps);
}

inline 
int TimeIntegrator::getIntegratorStep() const
{
   return(d_integrator_step);
}

inline 
int TimeIntegrator::getMaxIntegratorSteps() const
{
   return(d_max_integrator_steps);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif
}
}

