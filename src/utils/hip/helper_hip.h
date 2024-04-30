#ifndef AMP_HipHelpers
#define AMP_HipHelpers

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "AMP/utils/UtilityMacros.h"
#include "AMP/utils/hip/helper_string.h"

#include "StackTrace/source_location.h"

#include <hip/hip_runtime.h>


#ifndef EXIT_WAIVED
    #define EXIT_WAIVED 2
#endif


#ifdef __DRIVER_TYPES_H__
    #ifndef DEVICE_RESET
        #define DEVICE_RESET hipDeviceReset();
    #endif
#else
    #ifndef DEVICE_RESET
        #define DEVICE_RESET
    #endif
#endif


namespace AMP::Utilities {
enum class MemoryType : uint8_t;
}

// Get the pointer type from hip
AMP::Utilities::MemoryType getHipMemoryType( const void *ptr );


// Get the name of a return code
template<typename T>
const char *hipGetName( T result );

// Check the return code
template<typename T>
void checkHipErrors( T result,
                      const StackTrace::source_location &source = SOURCE_LOCATION_CURRENT() );

// Get the last hip error
void getLastHipError( const char *errorMessage,
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

#ifdef __HIP_RUNTIME_H__
// General GPU Device HIP Initialization
int gpuDeviceInit( int devID );

// This function returns the best GPU (with maximum GFLOPS)
int gpuGetMaxGflopsDeviceId();


// Initialization code to find the best HIP Device
int findHipDevice( int argc, const char **argv );

// General check for HIP GPU SM Capabilities
bool checkHipCapabilities( int major_version, int minor_version );


#endif

// end of HIP Helper Functions


#endif
