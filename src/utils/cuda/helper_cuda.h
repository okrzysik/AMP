#ifndef AMP_CudaHelpers
#define AMP_CudaHelpers

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AMP/utils/UtilityMacros.h"
#include "AMP/utils/cuda/helper_string.h"

#include "StackTrace/source_location.h"

#include <cuda.h>
#include <cuda_runtime.h>

#ifndef EXIT_WAIVED
    #define EXIT_WAIVED 2
#endif


#ifdef __DRIVER_TYPES_H__
    #ifndef DEVICE_RESET
        #define DEVICE_RESET cudaDeviceReset();
    #endif
#else
    #ifndef DEVICE_RESET
        #define DEVICE_RESET
    #endif
#endif


// Get the name of a return code
template<typename T>
const char *cudaGetName( T result );

// Check the return code
template<typename T>
void checkCudaErrors( T result,
                      const StackTrace::source_location &source = SOURCE_LOCATION_CURRENT() );

// Get the last cuda error
void getLastCudaError( const char *errorMessage,
                       const StackTrace::source_location &source = SOURCE_LOCATION_CURRENT() );

#ifndef MAX
    #define MAX( a, b ) ( a > b ? a : b )
#endif

// Float To Int conversion
inline int ftoi( float value )
{
    return ( value >= 0 ? (int) ( value + 0.5 ) : (int) ( value - 0.5 ) );
}

// Beginning of GPU Architecture definitions
inline int _ConvertSMVer2Cores( int major, int minor )
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = { { 0x20, 32 },  // Fermi Generation (SM 2.0) GF100 class
                                        { 0x21, 48 },  // Fermi Generation (SM 2.1) GF10x class
                                        { 0x30, 192 }, // Kepler Generation (SM 3.0) GK10x class
                                        { 0x32, 192 }, // Kepler Generation (SM 3.2) GK10x class
                                        { 0x35, 192 }, // Kepler Generation (SM 3.5) GK11x class
                                        { 0x37, 192 }, // Kepler Generation (SM 3.7) GK21x class
                                        { 0x50, 128 }, // Maxwell Generation (SM 5.0) GM10x class
                                        { 0x52, 128 }, // Maxwell Generation (SM 5.2) GM20x class
                                        { -1, -1 } };

    int index = 0;

    while ( nGpuArchCoresPerSM[index].SM != -1 ) {
        if ( nGpuArchCoresPerSM[index].SM == ( ( major << 4 ) + minor ) ) {
            return nGpuArchCoresPerSM[index].Cores;
        }

        index++;
    }

    // If we don't find the values, we default use the previous one to run properly
    printf( "MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n",
            major,
            minor,
            nGpuArchCoresPerSM[index - 1].Cores );
    return nGpuArchCoresPerSM[index - 1].Cores;
}
// end of GPU Architecture definitions

#ifdef __CUDA_RUNTIME_H__
// General GPU Device CUDA Initialization
inline int gpuDeviceInit( int devID )
{
    int device_count;
    checkCudaErrors( cudaGetDeviceCount( &device_count ) );

    if ( device_count == 0 ) {
        fprintf( stderr, "gpuDeviceInit() CUDA error: no devices supporting CUDA.\n" );
        exit( EXIT_FAILURE );
    }

    if ( devID < 0 ) {
        devID = 0;
    }

    if ( devID > device_count - 1 ) {
        fprintf( stderr, "\n" );
        fprintf( stderr, ">> %d CUDA capable GPU device(s) detected. <<\n", device_count );
        fprintf( stderr, ">> gpuDeviceInit (-device=%d) is not a valid GPU device. <<\n", devID );
        fprintf( stderr, "\n" );
        return -devID;
    }

    cudaDeviceProp deviceProp;
    checkCudaErrors( cudaGetDeviceProperties( &deviceProp, devID ) );

    if ( deviceProp.computeMode == cudaComputeModeProhibited ) {
        fprintf( stderr,
                 "Error: device is running in <Compute Mode Prohibited>, no threads can "
                 "use ::cudaSetDevice().\n" );
        return -1;
    }

    if ( deviceProp.major < 1 ) {
        fprintf( stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n" );
        exit( EXIT_FAILURE );
    }

    checkCudaErrors( cudaSetDevice( devID ) );
    printf( "gpuDeviceInit() CUDA Device [%d]: \"%s\n", devID, deviceProp.name );

    return devID;
}

// This function returns the best GPU (with maximum GFLOPS)
inline int gpuGetMaxGflopsDeviceId()
{
    int current_device = 0, sm_per_multiproc = 0;
    int max_perf_device = 0;
    int device_count = 0, best_SM_arch = 0;
    int devices_prohibited = 0;

    unsigned long long max_compute_perf = 0;
    cudaDeviceProp deviceProp;
    cudaGetDeviceCount( &device_count );

    checkCudaErrors( cudaGetDeviceCount( &device_count ) );

    if ( device_count == 0 ) {
        fprintf( stderr, "gpuGetMaxGflopsDeviceId() CUDA error: no devices supporting CUDA.\n" );
        exit( EXIT_FAILURE );
    }

    // Find the best major SM Architecture GPU device
    while ( current_device < device_count ) {
        cudaGetDeviceProperties( &deviceProp, current_device );

        // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
        if ( deviceProp.computeMode != cudaComputeModeProhibited ) {
            if ( deviceProp.major > 0 && deviceProp.major < 9999 ) {
                best_SM_arch = MAX( best_SM_arch, deviceProp.major );
            }
        } else {
            devices_prohibited++;
        }

        current_device++;
    }

    if ( devices_prohibited == device_count ) {
        fprintf(
            stderr,
            "gpuGetMaxGflopsDeviceId() CUDA error: all devices have compute mode prohibited.\n" );
        exit( EXIT_FAILURE );
    }

    // Find the best CUDA capable GPU device
    current_device = 0;

    while ( current_device < device_count ) {
        cudaGetDeviceProperties( &deviceProp, current_device );

        // If this GPU is not running on Compute Mode prohibited, then we can add it to the list
        if ( deviceProp.computeMode != cudaComputeModeProhibited ) {
            if ( deviceProp.major == 9999 && deviceProp.minor == 9999 ) {
                sm_per_multiproc = 1;
            } else {
                sm_per_multiproc = _ConvertSMVer2Cores( deviceProp.major, deviceProp.minor );
            }

            unsigned long long compute_perf = (unsigned long long) deviceProp.multiProcessorCount *
                                              sm_per_multiproc * deviceProp.clockRate;

            if ( compute_perf > max_compute_perf ) {
                // If we find GPU with SM major > 2, search only these
                if ( best_SM_arch > 2 ) {
                    // If our device==dest_SM_arch, choose this, or else pass
                    if ( deviceProp.major == best_SM_arch ) {
                        max_compute_perf = compute_perf;
                        max_perf_device  = current_device;
                    }
                } else {
                    max_compute_perf = compute_perf;
                    max_perf_device  = current_device;
                }
            }
        }

        ++current_device;
    }
    return max_perf_device;
}


// Initialization code to find the best CUDA Device
inline int findCudaDevice( int argc, const char **argv )
{
    cudaDeviceProp deviceProp;
    int devID = 0;

    // If the command-line has a device number specified, use it
    if ( checkCmdLineFlag( argc, argv, "device" ) ) {
        devID = getCmdLineArgumentInt( argc, argv, "device=" );

        if ( devID < 0 ) {
            printf( "Invalid command line parameter\n " );
            exit( EXIT_FAILURE );
        } else {
            devID = gpuDeviceInit( devID );

            if ( devID < 0 ) {
                printf( "exiting...\n" );
                exit( EXIT_FAILURE );
            }
        }
    } else {
        // Otherwise pick the device with highest Gflops/s
        devID = gpuGetMaxGflopsDeviceId();
        checkCudaErrors( cudaSetDevice( devID ) );
        checkCudaErrors( cudaGetDeviceProperties( &deviceProp, devID ) );
        printf( "GPU Device %d: \"%s\" with compute capability %d.%d\n\n",
                devID,
                deviceProp.name,
                deviceProp.major,
                deviceProp.minor );
    }

    return devID;
}

// General check for CUDA GPU SM Capabilities
inline bool checkCudaCapabilities( int major_version, int minor_version )
{
    cudaDeviceProp deviceProp;
    deviceProp.major = 0;
    deviceProp.minor = 0;
    int dev;

    checkCudaErrors( cudaGetDevice( &dev ) );
    checkCudaErrors( cudaGetDeviceProperties( &deviceProp, dev ) );

    if ( ( deviceProp.major > major_version ) ||
         ( deviceProp.major == major_version && deviceProp.minor >= minor_version ) ) {
        printf( "  Device %d: <%16s >, Compute SM %d.%d detected\n",
                dev,
                deviceProp.name,
                deviceProp.major,
                deviceProp.minor );
        return true;
    } else {
        printf( "  No GPU device was found that can support CUDA compute capability %d.%d.\n",
                major_version,
                minor_version );
        return false;
    }
}
#endif

// end of CUDA Helper Functions


#endif
