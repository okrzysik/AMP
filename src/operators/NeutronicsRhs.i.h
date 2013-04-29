#ifndef operator_NeutronicsRhs_i_h
#define operator_NeutronicsRhs_i_h

namespace AMP {
namespace Operator {


inline void NeutronicsRhs::setTimeInSeconds( double setSeconds ) {
    AMP_ASSERT (setSeconds >= 0.);
    // first assume time does not go backwards.
    int    timeStep = d_timeStep;
    double seconds  = d_timeStepInSeconds;
    if (setSeconds >= d_timeStepInSeconds ) {
      for( int i=d_timeStep; i < d_numTimeSteps; i++ ){ 
        seconds += d_timeStepsInDays[i] * d_secondsPerDay;
        if ( setSeconds > seconds ) {timeStep = i+1; }
        else {
          AMP_ASSERT(timeStep < d_numTimeSteps);
          d_timeStep = timeStep;
          d_timeStepInSeconds = seconds - d_timeStepsInDays[d_timeStep];
          return;
        }
      }
    } else {
      seconds  = 0.;
      timeStep = 0;
      for( int i=0; i < d_numTimeSteps; i++ ){ 
        seconds += d_timeStepsInDays[i] * d_secondsPerDay;
        if ( setSeconds > seconds ) {timeStep = i+1; }
        else {
          AMP_ASSERT(timeStep < d_numTimeSteps);
          d_timeStep = timeStep;
          d_timeStepInSeconds = seconds - d_timeStepsInDays[d_timeStep] * d_secondsPerDay;
          return;
        }
      }
    }
    AMP_INSIST(false,"Could not find the appropriate time.");
}
inline void NeutronicsRhs::setTimeInDays( double days ) {
    setTimeInSeconds(days*d_secondsPerDay);
}


}
} // end namespace AMP

#endif // operator_NeutronicsRhs_i_h
