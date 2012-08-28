#ifndef included_AMP_PetscMonitor
#define included_AMP_PetscMonitor

#include <set>
#include <map>
#include <complex>
#include "Utilities.h"
#include "AMP_MPI.h"

#ifdef USE_EXT_PETSC
    #include "petsc.h"
    #include "petscksp.h"
    #include "petscsnes.h"
#endif


namespace AMP {


/**
 * \class PetscMonitor
 *
 * @brief Provide a monitor that petsc can use instead of the default.
 * @details  This class provides routines for convergenece monitoring.
 *   The default convergence monitoring by petsc does not properly cleanup
 *   until PetscFinialize and creates problems when used on sub-communicators.
 */
class PetscMonitor
{
public:


    /**
     *\brief  Default constructor
     *\details  This is the default constructor that should be used
     *\param comm  Communicator to use for the monitor
     */
    PetscMonitor( AMP::AMP_MPI comm );

    //!  Empty deconstructor
    ~PetscMonitor();

    #ifdef USE_EXT_PETSC
        //! Routine to pass to petsc for monitoring KSP
        static PetscErrorCode monitorKSP(KSP,int,double,void*);

        //! Routine to pass to petsc for monitoring SNES
        static PetscErrorCode monitorSNES(SNES,int,double,void*);

        //! Routine to pass to petsc for monitoring KSP delete
        //static PetscErrorCode (*)(void**)  getKSPMonitorDelete();

        //! Routine to pass to petsc for monitoring SNES delete
        //static PetscErrorCode (*)(void**)  getSNESMonitorDelete();
    #endif

private:

    PetscMonitor();

    AMP::AMP_MPI d_comm;


};


}

#endif
