
#ifdef DEBUG_NO_INLINE
#define inline
#endif
namespace AMP{
namespace TimeIntegrator{

template<int DIM> inline
double ImplicitTimeIntegrator<DIM>::getCurrentDt() const
{
   return(d_current_dt);
}

template<int DIM> inline
double ImplicitTimeIntegrator<DIM>::getInitialTime() const
{
   return(d_initial_time);
}

template<int DIM> inline
double ImplicitTimeIntegrator<DIM>::getFinalTime() const
{
   return(d_final_time);
}

template<int DIM> inline
double ImplicitTimeIntegrator<DIM>::getCurrentTime() const
{
   return (d_current_time);
}

template<int DIM> inline
bool ImplicitTimeIntegrator<DIM>::stepsRemaining() const
{
   return(d_integrator_step < d_max_integrator_steps);
}

template<int DIM> inline 
int ImplicitTimeIntegrator<DIM>::getIntegratorStep() const
{
   return(d_integrator_step);
}

template<int DIM> inline 
int ImplicitTimeIntegrator<DIM>::getMaxIntegratorSteps() const
{
   return(d_max_integrator_steps);
}

#ifdef DEBUG_NO_INLINE
#undef inline
#endif

}
}

