# FindAMP
# ---------
#
# Find the Advanced Multi-Physics package (AMP)
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



# Set CMake policies
CMAKE_POLICY( SET CMP0057 NEW )


# Check that PROJ and ${PROJ}_INSTALL_DIR have been set
IF ( NOT PROJ )
    MESSAGE( FATAL_ERROR "PROJ must be set before calling FindTPLs")
ENDIF()
IF ( NOT ${PROJ}_INSTALL_DIR )
    MESSAGE( FATAL_ERROR "${PROJ}_INSTALL_DIR must be set before calling FindTPLs")
ENDIF()


# Disable link-time optimization
SET( DISABLE_LTO @DISABLE_LTO@ )


# Get the list of TPLs from AMP
SET( TPLS_REQUIRED @TPL_LIST_FOUND@ )
SET( TPLS_OPTIONAL )


# Add any addional TPLs needed
FOREACH( tmp ${EXTRA_REQUIRED_TPLS} )
    IF ( NOT ( ${tmp} IN_LIST TPLS_REQUIRED ) )
        SET( TPLS_REQUIRED ${TPLS_REQUIRED} ${tmp} )
    ENDIF()
ENDFOREACH()
FOREACH( tmp ${EXTRA_OPTIONAL_TPLS} )
    IF ( NOT ( ${tmp} IN_LIST TPLS_REQUIRED ) AND NOT ( ${tmp} IN_LIST TPLS_OPTIONAL ) )
        SET( TPLS_OPTIONAL ${TPLS_OPTIONAL} ${tmp} )
    ENDIF()
ENDFOREACH()


# Call the TPL builder with the appropriate TPLs for the AMP
SET( TPL_DIRECTORY "@TPL_DIRECTORY@" )
SET( CMAKE_MODULE_PATH "@TPL_DIRECTORY@" ${CMAKE_MODULE_PATH} )
FIND_PACKAGE( TPLs REQUIRED ${TPLS_REQUIRED} OPTIONAL_COMPONENTS ${TPLS_OPTIONAL} )
FOREACH( tmp ${TPL_LIST} )
    IF ( TPL_FOUND_${tmp} )
        SET( TPL_LIST_FOUND ${TPL_LIST_FOUND} ${tmp} )
    ENDIF()
ENDFOREACH()


# Load a dummy timer if one was not include
IF ( NOT DEFINED TIMER_INCLUDE )
    INCLUDE( "${TPL_DIRECTORY}/cmake/FindTimer.cmake" )
    CONFIGURE_TIMER( FALSE "${${PROJ}_INSTALL_DIR}/include" TRUE )
    SET( TPL_INCLUDE_DIRS ${TPL_INCLUDE_DIRS} ${TIMER_INCLUDE} )
ENDIF()


# Set the maximum number of processors for a test
IF ( NOT TEST_MAX_PROCS )
    SET( TEST_MAX_PROCS @TEST_MAX_PROCS@ )
ENDIF()

# Override the compiler flags (to add extra flags set by AMP)
SET( CMAKE_C_FLAGS "@CMAKE_C_FLAGS@" )
SET( CMAKE_CXX_FLAGS "@CMAKE_CXX_FLAGS@" )
SET( CMAKE_Fortran_FLAGS "@CMAKE_Fortran_FLAGS@" )
SET( LDLIBS "@LDLIBS@" )
SET( LDFLAGS "@LDFLAGS@" )
SET( LDLIBS_EXTRA "@LDLIBS_EXTRA@" )
SET( LDFLAGS_EXTRA "@LDFLAGS_EXTRA@" )


# Add the libraries for AMP
SET( AMP_FOUND TRUE )
SET( CMAKE_INSTALL_RPATH @CMAKE_INSTALL_RPATH@ ${CMAKE_INSTALL_RPATH} )
FIND_LIBRARY( AMP_LIB  NAMES @AMP_LIB@  PATHS "@AMP_INSTALL_DIR@/lib" NO_DEFAULT_PATH )
SET( AMP_LIBRARIES ${AMP_LIB} )
SET( AMP_INCLUDE_DIRS "@AMP_INSTALL_DIR@/include" )
SET( AMP_SOURCE_DIR "@AMP_SOURCE_DIR@" )
SET( AMP_MACRO_CMAKE "${TPL_MACRO_CMAKE}" )


# Setup doxygen
SET( USE_DOXYGEN @USE_DOXYGEN@ )
SET( USE_EXT_DOXYGEN @USE_EXT_DOXYGEN@ )
IF ( USE_DOXYGEN OR USE_EXT_DOXYGEN )
    SET( USE_DOXYGEN TRUE )
    SET( USE_EXT_DOXYGEN TRUE )
ENDIF()


# Add documentation folders
SET( AMP_DOC_DIRS @AMP_DOC_DIRS@ )



