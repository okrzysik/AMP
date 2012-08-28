#include "PetscMonitor.h"

namespace AMP {


/********************************************************************
*  Constructor                                                      *
********************************************************************/
PetscMonitor::PetscMonitor( AMP_MPI comm ) 
{
    d_comm = comm;
}


/********************************************************************
*  Routines to provide petsc with function pointers for monitoring  *
********************************************************************/
#ifdef USE_EXT_PETSC
PetscErrorCode PetscMonitor::monitorKSP(KSP,int,double,void*)
{
}
PetscErrorCode PetscMonitor::monitorSNES(SNES,int,double,void*)
{
}
#endif


}

