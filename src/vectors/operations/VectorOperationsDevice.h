#ifndef included_AMP_VectorOperationsDevice
#define included_AMP_VectorOperationsDevice

#ifdef USE_CUDA
    #include "AMP/vectors/operations/cuda/VectorOperationsCuda.h"
#endif
#ifdef USE_HIP
    #include "AMP/vectors/operations/hip/VectorOperationsHip.h"
#endif

#endif
