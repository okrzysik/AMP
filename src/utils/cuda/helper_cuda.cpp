#include "AMP/utils/cuda/helper_cuda.h"

#include <cuda.h>
#include <cuda_runtime.h>


// Check
template<typename T>
void check( T result, char const *const func, const char *const file, int const line )
{
<<<<<<< HEAD
    switch ( error ) {
    case cudaSuccess:
        return "cudaSuccess";
    case cudaErrorMissingConfiguration:
        return "cudaErrorMissingConfiguration";
    case cudaErrorMemoryAllocation:
        return "cudaErrorMemoryAllocation";
    case cudaErrorInitializationError:
        return "cudaErrorInitializationError";
    case cudaErrorLaunchFailure:
        return "cudaErrorLaunchFailure";
    case cudaErrorPriorLaunchFailure:
        return "cudaErrorPriorLaunchFailure";
    case cudaErrorLaunchTimeout:
        return "cudaErrorLaunchTimeout";
    case cudaErrorLaunchOutOfResources:
        return "cudaErrorLaunchOutOfResources";
    case cudaErrorInvalidDeviceFunction:
        return "cudaErrorInvalidDeviceFunction";
    case cudaErrorInvalidConfiguration:
        return "cudaErrorInvalidConfiguration";
    case cudaErrorInvalidDevice:
        return "cudaErrorInvalidDevice";
    case cudaErrorInvalidValue:
        return "cudaErrorInvalidValue";
    case cudaErrorInvalidPitchValue:
        return "cudaErrorInvalidPitchValue";
    case cudaErrorInvalidSymbol:
        return "cudaErrorInvalidSymbol";
    case cudaErrorMapBufferObjectFailed:
        return "cudaErrorMapBufferObjectFailed";
    case cudaErrorUnmapBufferObjectFailed:
        return "cudaErrorUnmapBufferObjectFailed";
    case cudaErrorInvalidHostPointer:
        return "cudaErrorInvalidHostPointer";
    case cudaErrorInvalidDevicePointer:
        return "cudaErrorInvalidDevicePointer";
    case cudaErrorInvalidTexture:
        return "cudaErrorInvalidTexture";
    case cudaErrorInvalidTextureBinding:
        return "cudaErrorInvalidTextureBinding";
    case cudaErrorInvalidChannelDescriptor:
        return "cudaErrorInvalidChannelDescriptor";
    case cudaErrorInvalidMemcpyDirection:
        return "cudaErrorInvalidMemcpyDirection";
    case cudaErrorAddressOfConstant:
        return "cudaErrorAddressOfConstant";
    case cudaErrorTextureFetchFailed:
        return "cudaErrorTextureFetchFailed";
    case cudaErrorTextureNotBound:
        return "cudaErrorTextureNotBound";
    case cudaErrorSynchronizationError:
        return "cudaErrorSynchronizationError";
    case cudaErrorInvalidFilterSetting:
        return "cudaErrorInvalidFilterSetting";
    case cudaErrorInvalidNormSetting:
        return "cudaErrorInvalidNormSetting";
    case cudaErrorMixedDeviceExecution:
        return "cudaErrorMixedDeviceExecution";
    case cudaErrorCudartUnloading:
        return "cudaErrorCudartUnloading";
    case cudaErrorUnknown:
        return "cudaErrorUnknown";
    case cudaErrorNotYetImplemented:
        return "cudaErrorNotYetImplemented";
    case cudaErrorMemoryValueTooLarge:
        return "cudaErrorMemoryValueTooLarge";
    case cudaErrorInvalidResourceHandle:
        return "cudaErrorInvalidResourceHandle";
    case cudaErrorNotReady:
        return "cudaErrorNotReady";
    case cudaErrorInsufficientDriver:
        return "cudaErrorInsufficientDriver";
    case cudaErrorSetOnActiveProcess:
        return "cudaErrorSetOnActiveProcess";
    case cudaErrorInvalidSurface:
        return "cudaErrorInvalidSurface";
    case cudaErrorNoDevice:
        return "cudaErrorNoDevice";
    case cudaErrorECCUncorrectable:
        return "cudaErrorECCUncorrectable";
    case cudaErrorSharedObjectSymbolNotFound:
        return "cudaErrorSharedObjectSymbolNotFound";
    case cudaErrorSharedObjectInitFailed:
        return "cudaErrorSharedObjectInitFailed";
    case cudaErrorUnsupportedLimit:
        return "cudaErrorUnsupportedLimit";
    case cudaErrorDuplicateVariableName:
        return "cudaErrorDuplicateVariableName";
    case cudaErrorDuplicateTextureName:
        return "cudaErrorDuplicateTextureName";
    case cudaErrorDuplicateSurfaceName:
        return "cudaErrorDuplicateSurfaceName";
    case cudaErrorDevicesUnavailable:
        return "cudaErrorDevicesUnavailable";
    case cudaErrorInvalidKernelImage:
        return "cudaErrorInvalidKernelImage";
    case cudaErrorNoKernelImageForDevice:
        return "cudaErrorNoKernelImageForDevice";
    case cudaErrorIncompatibleDriverContext:
        return "cudaErrorIncompatibleDriverContext";
    case cudaErrorPeerAccessAlreadyEnabled:
        return "cudaErrorPeerAccessAlreadyEnabled";
    case cudaErrorPeerAccessNotEnabled:
        return "cudaErrorPeerAccessNotEnabled";
    case cudaErrorDeviceAlreadyInUse:
        return "cudaErrorDeviceAlreadyInUse";
    case cudaErrorProfilerDisabled:
        return "cudaErrorProfilerDisabled";
    case cudaErrorProfilerNotInitialized:
        return "cudaErrorProfilerNotInitialized";
    case cudaErrorProfilerAlreadyStarted:
        return "cudaErrorProfilerAlreadyStarted";
    case cudaErrorProfilerAlreadyStopped:
        return "cudaErrorProfilerAlreadyStopped";
    case cudaErrorAssert:
        return "cudaErrorAssert";
    case cudaErrorTooManyPeers:
        return "cudaErrorTooManyPeers";
    case cudaErrorHostMemoryAlreadyRegistered:
        return "cudaErrorHostMemoryAlreadyRegistered";
    case cudaErrorHostMemoryNotRegistered:
        return "cudaErrorHostMemoryNotRegistered";
    case cudaErrorOperatingSystem:
        return "cudaErrorOperatingSystem";
    case cudaErrorPeerAccessUnsupported:
        return "cudaErrorPeerAccessUnsupported";
    case cudaErrorLaunchMaxDepthExceeded:
        return "cudaErrorLaunchMaxDepthExceeded";
    case cudaErrorLaunchFileScopedTex:
        return "cudaErrorLaunchFileScopedTex";
    case cudaErrorLaunchFileScopedSurf:
        return "cudaErrorLaunchFileScopedSurf";
    case cudaErrorSyncDepthExceeded:
        return "cudaErrorSyncDepthExceeded";
    case cudaErrorLaunchPendingCountExceeded:
        return "cudaErrorLaunchPendingCountExceeded";
    case cudaErrorNotPermitted:
        return "cudaErrorNotPermitted";
    case cudaErrorNotSupported:
        return "cudaErrorNotSupported";
    case cudaErrorHardwareStackError:
        return "cudaErrorHardwareStackError";
    case cudaErrorIllegalInstruction:
        return "cudaErrorIllegalInstruction";
    case cudaErrorMisalignedAddress:
        return "cudaErrorMisalignedAddress";
    case cudaErrorInvalidAddressSpace:
        return "cudaErrorInvalidAddressSpace";
    case cudaErrorInvalidPc:
        return "cudaErrorInvalidPc";
    case cudaErrorIllegalAddress:
        return "cudaErrorIllegalAddress";
    case cudaErrorInvalidPtx:
        return "cudaErrorInvalidPtx";
    case cudaErrorInvalidGraphicsContext:
        return "cudaErrorInvalidGraphicsContext";
    case cudaErrorStartupFailure:
        return "cudaErrorStartupFailure";
    case cudaErrorApiFailureBase:
        return "cudaErrorApiFailureBase";
    case cudaErrorNvlinkUncorrectable:
        return "cudaErrorNvlinkUncorrectable";
    case cudaErrorJitCompilerNotFound:
        return "cudaErrorJitCompilerNotFound";
    case cudaErrorCooperativeLaunchTooLarge:
        return "cudaErrorCooperativeLaunchTooLarge";
    case cudaErrorDeviceUninitialized:
        return "cudaErrorDeviceUninitialized";
    case cudaErrorArrayIsMapped:
        return "cudaErrorArrayIsMapped";
    case cudaErrorAlreadyMapped:
        return "cudaErrorAlreadyMapped";
    case cudaErrorAlreadyAcquired:
        return "cudaErrorAlreadyAcquired";
    case cudaErrorNotMapped:
        return "cudaErrorNotMapped";
    case cudaErrorNotMappedAsArray:
        return "cudaErrorNotMappedAsArray";
    case cudaErrorNotMappedAsPointer:
        return "cudaErrorNotMappedAsPointer";
    case cudaErrorInvalidSource:
        return "cudaErrorInvalidSource";
    case cudaErrorFileNotFound:
        return "cudaErrorFileNotFound";
    case cudaErrorIllegalState:
        return "cudaErrorIllegalState";
    case cudaErrorSymbolNotFound:
        return "cudaErrorSymbolNotFound";
    case cudaErrorLaunchIncompatibleTexturing:
        return "cudaErrorLaunchIncompatibleTexturing";
    case cudaErrorContextIsDestroyed:
        return "cudaErrorContextIsDestroyed";
    case cudaErrorSystemNotReady:
        return "cudaErrorSystemNotReady";
    case cudaErrorSystemDriverMismatch:
        return "cudaErrorSystemDriverMismatch";
    case cudaErrorCompatNotSupportedOnDevice:
        return "cudaErrorCompatNotSupportedOnDevice";
    case cudaErrorStreamCaptureUnsupported:
        return "cudaErrorStreamCaptureUnsupported";
    case cudaErrorStreamCaptureInvalidated:
        return "cudaErrorStreamCaptureInvalidated";
    case cudaErrorStreamCaptureMerge:
        return "cudaErrorStreamCaptureMerge";
    case cudaErrorStreamCaptureUnmatched:
        return "cudaErrorStreamCaptureUnmatched";
    case cudaErrorStreamCaptureUnjoined:
        return "cudaErrorStreamCaptureUnjoined";
    case cudaErrorStreamCaptureIsolation:
        return "cudaErrorStreamCaptureIsolation";
    case cudaErrorStreamCaptureImplicit:
        return "cudaErrorStreamCaptureImplicit";
    case cudaErrorCapturedEvent:
        return "cudaErrorCapturedEvent";
    case cudaErrorStreamCaptureWrongThread:
        return "cudaErrorStreamCaptureWrongThread";
=======
    if ( result ) {
        fprintf( stderr,
                 "CUDA error at %s:%d code=%d(%s) \"%s\" \n",
                 file,
                 line,
                 static_cast<int>( result ),
                 cudaGetName( result ),
                 func );
        // Make sure we call CUDA Device Reset before exiting
        DEVICE_RESET
        exit( EXIT_FAILURE );
>>>>>>> 4c89191beaa0059acd51b51869ace29bbd4f0d32
    }
}


// CUDA Runtime error messages
template<>
const char *cudaGetName<cudaError_t>( cudaError_t error )
{
    return cudaGetErrorName( error );
}
template void check<cudaError_t>( cudaError_t, char const *const, const char *const, int const );


// CUDA Driver API errors
template<>
const char *cudaGetName<CUresult>( CUresult error )
{
    const char *str = nullptr;
    cuGetErrorName( error, &str );
    return str;
}
template void check<CUresult>( CUresult, char const *const, const char *const, int const );


// cuBLAS API errors
#ifdef CUBLAS_API_H_
template<>
const char *cudaGetName<cublasStatus_t>( cublasStatus_t error )
{
    switch ( error ) {
    case CUBLAS_STATUS_SUCCESS:
        return "CUBLAS_STATUS_SUCCESS";
    case CUBLAS_STATUS_NOT_INITIALIZED:
        return "CUBLAS_STATUS_NOT_INITIALIZED";
    case CUBLAS_STATUS_ALLOC_FAILED:
        return "CUBLAS_STATUS_ALLOC_FAILED";
    case CUBLAS_STATUS_INVALID_VALUE:
        return "CUBLAS_STATUS_INVALID_VALUE";
    case CUBLAS_STATUS_ARCH_MISMATCH:
        return "CUBLAS_STATUS_ARCH_MISMATCH";
    case CUBLAS_STATUS_MAPPING_ERROR:
        return "CUBLAS_STATUS_MAPPING_ERROR";
    case CUBLAS_STATUS_EXECUTION_FAILED:
        return "CUBLAS_STATUS_EXECUTION_FAILED";
    case CUBLAS_STATUS_INTERNAL_ERROR:
        return "CUBLAS_STATUS_INTERNAL_ERROR";
    case CUBLAS_STATUS_NOT_SUPPORTED:
        return "CUBLAS_STATUS_NOT_SUPPORTED";
    case CUBLAS_STATUS_LICENSE_ERROR:
        return "CUBLAS_STATUS_LICENSE_ERROR";
    }
    return "<unknown>";
}
template void
check<cublasStatus_t>( cublasStatus_t, char const *const, const char *const, int const );
#endif


// cuFFT API errors
#ifdef _CUFFT_H_
template<>
const char *cudaGetName<cufftResult>( cufftResult error )
{
    switch ( error ) {
    case CUFFT_SUCCESS:
        return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN:
        return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED:
        return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE:
        return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE:
        return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR:
        return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED:
        return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED:
        return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE:
        return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA:
        return "CUFFT_UNALIGNED_DATA";
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
        return "CUFFT_INCOMPLETE_PARAMETER_LIST";
    case CUFFT_INVALID_DEVICE:
        return "CUFFT_INVALID_DEVICE";
    case CUFFT_PARSE_ERROR:
        return "CUFFT_PARSE_ERROR";
    case CUFFT_NO_WORKSPACE:
        return "CUFFT_NO_WORKSPACE";
    case CUFFT_NOT_IMPLEMENTED:
        return "CUFFT_NOT_IMPLEMENTED";
    case CUFFT_LICENSE_ERROR:
        return "CUFFT_LICENSE_ERROR";
    }
    return "<unknown>";
}
template void check<cufftResult>( cufftResult, char const *const, const char *const, int const );
#endif


// cuSPARSE API errors
#ifdef CUSPARSEAPI
template<>
const char *cudaGetName<cusparseStatus_t>( cusparseStatus_t error )
{
    switch ( error ) {
    case CUSPARSE_STATUS_SUCCESS:
        return "CUSPARSE_STATUS_SUCCESS";
    case CUSPARSE_STATUS_NOT_INITIALIZED:
        return "CUSPARSE_STATUS_NOT_INITIALIZED";
    case CUSPARSE_STATUS_ALLOC_FAILED:
        return "CUSPARSE_STATUS_ALLOC_FAILED";
    case CUSPARSE_STATUS_INVALID_VALUE:
        return "CUSPARSE_STATUS_INVALID_VALUE";
    case CUSPARSE_STATUS_ARCH_MISMATCH:
        return "CUSPARSE_STATUS_ARCH_MISMATCH";
    case CUSPARSE_STATUS_MAPPING_ERROR:
        return "CUSPARSE_STATUS_MAPPING_ERROR";
    case CUSPARSE_STATUS_EXECUTION_FAILED:
        return "CUSPARSE_STATUS_EXECUTION_FAILED";
    case CUSPARSE_STATUS_INTERNAL_ERROR:
        return "CUSPARSE_STATUS_INTERNAL_ERROR";
    case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    }
    return "<unknown>";
}
template void
check<cusparseStatus_t>( cufftResult, char const *const, const char *const, int const );
#endif


// cuSOLVER API errors
#ifdef CUSOLVER_COMMON_H_
template<>
const char *cudaGetName<cusolverStatus_t>( cusolverStatus_t error )
{
    switch ( error ) {
    case CUSOLVER_STATUS_SUCCESS:
        return "CUSOLVER_STATUS_SUCCESS";
    case CUSOLVER_STATUS_NOT_INITIALIZED:
        return "CUSOLVER_STATUS_NOT_INITIALIZED";
    case CUSOLVER_STATUS_ALLOC_FAILED:
        return "CUSOLVER_STATUS_ALLOC_FAILED";
    case CUSOLVER_STATUS_INVALID_VALUE:
        return "CUSOLVER_STATUS_INVALID_VALUE";
    case CUSOLVER_STATUS_ARCH_MISMATCH:
        return "CUSOLVER_STATUS_ARCH_MISMATCH";
    case CUSOLVER_STATUS_MAPPING_ERROR:
        return "CUSOLVER_STATUS_MAPPING_ERROR";
    case CUSOLVER_STATUS_EXECUTION_FAILED:
        return "CUSOLVER_STATUS_EXECUTION_FAILED";
    case CUSOLVER_STATUS_INTERNAL_ERROR:
        return "CUSOLVER_STATUS_INTERNAL_ERROR";
    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    case CUSOLVER_STATUS_NOT_SUPPORTED:
        return "CUSOLVER_STATUS_NOT_SUPPORTED ";
    case CUSOLVER_STATUS_ZERO_PIVOT:
        return "CUSOLVER_STATUS_ZERO_PIVOT";
    case CUSOLVER_STATUS_INVALID_LICENSE:
        return "CUSOLVER_STATUS_INVALID_LICENSE";
    }

    return "<unknown>";
}
template void
check<cusolverStatus_t>( cufftResult, char const *const, const char *const, int const );
#endif


// cuRAND API errors
#ifdef CURAND_H_
template<>
const char *cudaGetName<curandStatus_t>( curandStatus_t error )
{
    switch ( error ) {
    case CURAND_STATUS_SUCCESS:
        return "CURAND_STATUS_SUCCESS";
    case CURAND_STATUS_VERSION_MISMATCH:
        return "CURAND_STATUS_VERSION_MISMATCH";
    case CURAND_STATUS_NOT_INITIALIZED:
        return "CURAND_STATUS_NOT_INITIALIZED";
    case CURAND_STATUS_ALLOCATION_FAILED:
        return "CURAND_STATUS_ALLOCATION_FAILED";
    case CURAND_STATUS_TYPE_ERROR:
        return "CURAND_STATUS_TYPE_ERROR";
    case CURAND_STATUS_OUT_OF_RANGE:
        return "CURAND_STATUS_OUT_OF_RANGE";
    case CURAND_STATUS_LENGTH_NOT_MULTIPLE:
        return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
    case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED:
        return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
    case CURAND_STATUS_LAUNCH_FAILURE:
        return "CURAND_STATUS_LAUNCH_FAILURE";
    case CURAND_STATUS_PREEXISTING_FAILURE:
        return "CURAND_STATUS_PREEXISTING_FAILURE";
    case CURAND_STATUS_INITIALIZATION_FAILED:
        return "CURAND_STATUS_INITIALIZATION_FAILED";
    case CURAND_STATUS_ARCH_MISMATCH:
        return "CURAND_STATUS_ARCH_MISMATCH";
    case CURAND_STATUS_INTERNAL_ERROR:
        return "CURAND_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
}
template void check<curandStatus_t>( cufftResult, char const *const, const char *const, int const );
#endif


// NPP API errors
#ifdef NV_NPPIDEFS_H
template<>
const char *cudaGetName<NppStatus>( NppStatus error )
{
    switch ( error ) {
    case NPP_NOT_SUPPORTED_MODE_ERROR:
        return "NPP_NOT_SUPPORTED_MODE_ERROR";
    case NPP_ROUND_MODE_NOT_SUPPORTED_ERROR:
        return "NPP_ROUND_MODE_NOT_SUPPORTED_ERROR";
    case NPP_RESIZE_NO_OPERATION_ERROR:
        return "NPP_RESIZE_NO_OPERATION_ERROR";
    case NPP_NOT_SUFFICIENT_COMPUTE_CAPABILITY:
        return "NPP_NOT_SUFFICIENT_COMPUTE_CAPABILITY";
    case NPP_BAD_ARGUMENT_ERROR:
        return "NPP_BAD_ARGUMENT_ERROR";
    case NPP_COEFFICIENT_ERROR:
        return "NPP_COEFFICIENT_ERROR";
    case NPP_RECTANGLE_ERROR:
        return "NPP_RECTANGLE_ERROR";
    case NPP_QUADRANGLE_ERROR:
        return "NPP_QUADRANGLE_ERROR";
    case NPP_MEMORY_ALLOCATION_ERR:
        return "NPP_MEMORY_ALLOCATION_ERROR";
    case NPP_HISTOGRAM_NUMBER_OF_LEVELS_ERROR:
        return "NPP_HISTOGRAM_NUMBER_OF_LEVELS_ERROR";
    case NPP_INVALID_HOST_POINTER_ERROR:
        return "NPP_INVALID_HOST_POINTER_ERROR";
    case NPP_INVALID_DEVICE_POINTER_ERROR:
        return "NPP_INVALID_DEVICE_POINTER_ERROR";
    case NPP_LUT_NUMBER_OF_LEVELS_ERROR:
        return "NPP_LUT_NUMBER_OF_LEVELS_ERROR";
    case NPP_TEXTURE_BIND_ERROR:
        return "NPP_TEXTURE_BIND_ERROR";
    case NPP_WRONG_INTERSECTION_ROI_ERROR:
        return "NPP_WRONG_INTERSECTION_ROI_ERROR";
    case NPP_NOT_EVEN_STEP_ERROR:
        return "NPP_NOT_EVEN_STEP_ERROR";
    case NPP_INTERPOLATION_ERROR:
        return "NPP_INTERPOLATION_ERROR";
    case NPP_RESIZE_FACTOR_ERROR:
        return "NPP_RESIZE_FACTOR_ERROR";
    case NPP_HAAR_CLASSIFIER_PIXEL_MATCH_ERROR:
        return "NPP_HAAR_CLASSIFIER_PIXEL_MATCH_ERROR";
    case NPP_MEMFREE_ERROR:
        return "NPP_MEMFREE_ERROR";
    case NPP_MEMSET_ERROR:
        return "NPP_MEMSET_ERROR";
    case NPP_MEMCPY_ERROR:
        return "NPP_MEMCPY_ERROR";
    case NPP_MIRROR_FLIP_ERROR:
        return "NPP_MIRROR_FLIP_ERROR";
    case NPP_ALIGNMENT_ERROR:
        return "NPP_ALIGNMENT_ERROR";
    case NPP_STEP_ERROR:
        return "NPP_STEP_ERROR";
    case NPP_SIZE_ERROR:
        return "NPP_SIZE_ERROR";
    case NPP_NULL_POINTER_ERROR:
        return "NPP_NULL_POINTER_ERROR";
    case NPP_CUDA_KERNEL_EXECUTION_ERROR:
        return "NPP_CUDA_KERNEL_EXECUTION_ERROR";
    case NPP_NOT_IMPLEMENTED_ERROR:
        return "NPP_NOT_IMPLEMENTED_ERROR";
    case NPP_ERROR:
        return "NPP_ERROR";
    case NPP_SUCCESS:
        return "NPP_SUCCESS";
    case NPP_WRONG_INTERSECTION_QUAD_WARNING:
        return "NPP_WRONG_INTERSECTION_QUAD_WARNING";
    case NPP_MISALIGNED_DST_ROI_WARNING:
        return "NPP_MISALIGNED_DST_ROI_WARNING";
    case NPP_AFFINE_QUAD_INCORRECT_WARNING:
        return "NPP_AFFINE_QUAD_INCORRECT_WARNING";
    case NPP_DOUBLE_SIZE_WARNING:
        return "NPP_DOUBLE_SIZE_WARNING";
    case NPP_WRONG_INTERSECTION_ROI_WARNING:
        return "NPP_WRONG_INTERSECTION_ROI_WARNING";
    case NPP_LUT_PALETTE_BITSIZE_ERROR:
        return "NPP_LUT_PALETTE_BITSIZE_ERROR";
    case NPP_ZC_MODE_NOT_SUPPORTED_ERROR:
        return "NPP_ZC_MODE_NOT_SUPPORTED_ERROR";
    case NPP_QUALITY_INDEX_ERROR:
        return "NPP_QUALITY_INDEX_ERROR";
    case NPP_CHANNEL_ORDER_ERROR:
        return "NPP_CHANNEL_ORDER_ERROR";
    case NPP_ZERO_MASK_VALUE_ERROR:
        return "NPP_ZERO_MASK_VALUE_ERROR";
    case NPP_NUMBER_OF_CHANNELS_ERROR:
        return "NPP_NUMBER_OF_CHANNELS_ERROR";
    case NPP_COI_ERROR:
        return "NPP_COI_ERROR";
    case NPP_DIVISOR_ERROR:
        return "NPP_DIVISOR_ERROR";
    case NPP_CHANNEL_ERROR:
        return "NPP_CHANNEL_ERROR";
    case NPP_STRIDE_ERROR:
        return "NPP_STRIDE_ERROR";
    case NPP_ANCHOR_ERROR:
        return "NPP_ANCHOR_ERROR";
    case NPP_MASK_SIZE_ERROR:
        return "NPP_MASK_SIZE_ERROR";
    case NPP_MOMENT_00_ZERO_ERROR:
        return "NPP_MOMENT_00_ZERO_ERROR";
    case NPP_THRESHOLD_NEGATIVE_LEVEL_ERROR:
        return "NPP_THRESHOLD_NEGATIVE_LEVEL_ERROR";
    case NPP_THRESHOLD_ERROR:
        return "NPP_THRESHOLD_ERROR";
    case NPP_CONTEXT_MATCH_ERROR:
        return "NPP_CONTEXT_MATCH_ERROR";
    case NPP_FFT_FLAG_ERROR:
        return "NPP_FFT_FLAG_ERROR";
    case NPP_FFT_ORDER_ERROR:
        return "NPP_FFT_ORDER_ERROR";
    case NPP_SCALE_RANGE_ERROR:
        return "NPP_SCALE_RANGE_ERROR";
    case NPP_DATA_TYPE_ERROR:
        return "NPP_DATA_TYPE_ERROR";
    case NPP_OUT_OFF_RANGE_ERROR:
        return "NPP_OUT_OFF_RANGE_ERROR";
    case NPP_DIVIDE_BY_ZERO_ERROR:
        return "NPP_DIVIDE_BY_ZERO_ERROR";
    case NPP_RANGE_ERROR:
        return "NPP_RANGE_ERROR";
    case NPP_NO_MEMORY_ERROR:
        return "NPP_NO_MEMORY_ERROR";
    case NPP_ERROR_RESERVED:
        return "NPP_ERROR_RESERVED";
    case NPP_NO_OPERATION_WARNING:
        cudaGetName return "NPP_NO_OPERATION_WARNING";
    case NPP_DIVIDE_BY_ZERO_WARNING:
        return "NPP_DIVIDE_BY_ZERO_WARNING";
    }
    return "<unknown>";
}
template void check<curandStatus_t>( NppStatus, char const *const, const char *const, int const );
#endif
