# FindAMP
# ---------
#
# Find the Advanved Multi-Physics package (AMP)
#
# Use this module by invoking find_package with the form:
#
#   find_package( AMP
#     [version] [EXACT]         # Minimum or EXACT version e.g. 1.36.0
#     [REQUIRED]                # Fail with error if the TPLs are not found
#     [COMPONENTS <libs>...]    # List of TPLs to include
#   )
#
# This module finds headers and requested component libraries for the AMP
#
#   AMP_FOUND           - True if headers and requested libraries were found
#   AMP_LIBRARIES       - AMP libraries (and dependencies)
#   AMP_INCLUDE_DIRS    - TPL Include paths
#   AMP_SOURCE_DIR      - Source directory   
#   AMP_MACRO_CMAKE     - File to macros.cmake provided by the TPL install

# Check that PROJ and ${PROJ}_INSTALL_DIR have been set
IF ( NOT PROJ )
    MESSAGE( FATAL_ERROR "PROJ must be set before calling FindTPLs")
ENDIF()
IF ( NOT ${PROJ}_INSTALL_DIR )
    MESSAGE( FATAL_ERROR "${PROJ}_INSTALL_DIR must be set before calling FindTPLs")
ENDIF()


# Call the TPL builder with the appropriate TPLs for the AMP
SET( CMAKE_MODULE_PATH "@TPL_DIRECTORY@" ${CMAKE_MODULE_PATH} )
FIND_PACKAGE( TPLs REQUIRED @TPL_LIST_FOUND@ OPTIONAL_COMPONENTS )


# Set the maximum number of processors for a test
IF ( NOT TEST_MAX_PROCS )
    SET( TEST_MAX_PROCS @TEST_MAX_PROCS@ )
ENDIF()


# Set MATLAB variables (eventually needs to be moved to the TPL builder)
SET( USE_MATLAB @USE_MATLAB@ )
SET( MATLAB_DIRECTORY @MATLAB_DIRECTORY@ )
SET( MEX_FLAGS @MEX_FLAGS@ )
SET( MEX_INCLUDE @MEX_INCLUDE@ )
SET( MEX_LDFLAGS @MEX_LDFLAGS@ )
SET( MEX_LIBS @MEX_LIBS@ )
SET( MATLAB_TARGET @MATLAB_TARGET@ )
SET( MEX_EXTENSION @MEX_EXTENSION@ )
SET( MATLAB_EXTERN @MATLAB_EXTERN@ )
SET( MATLAB_OS @MATLAB_OS@ )


# Add the libraries for AMP
SET( AMP_FOUND TRUE )
SET( CMAKE_INSTALL_RPATH @CMAKE_INSTALL_RPATH@ ${CMAKE_INSTALL_RPATH} )
ADD_DEFINITIONS( -DUSE_AMP_MODEL )
FIND_LIBRARY( AMP_LIB  NAMES @AMP_LIB@  PATHS "@AMP_INSTALL_DIR@/lib" NO_DEFAULT_PATH )
SET( AMP_LIBRARIES ${AMP_LIB} ${TPL_LIBRARIES} )
SET( AMP_INCLUDE_DIRS "@AMP_INSTALL_DIR@/include" ${TPL_INCLUDE_DIRS} )
SET( AMP_SOURCE_DIR "@AMP_SOURCE_DIR@" )
SET( AMP_MACRO_CMAKE "${TPL_MACRO_CMAKE}" )

