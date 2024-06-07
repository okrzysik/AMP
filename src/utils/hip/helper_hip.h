#ifndef AMP_HipHelpers
#define AMP_HipHelpers

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hip/hip_runtime.h>

#include "AMP/utils/UtilityMacros.h"
#include "StackTrace/source_location.h"

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

#endif
