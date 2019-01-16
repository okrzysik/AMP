# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 3.9)


# Print banner
MESSAGE ("===============")
MESSAGE ("Configuring AMP")
MESSAGE ("===============")


# Set the project name
SET( PROJ AMP )                 # Set the project name for CMake
SET( AMP_LIB amp )              # Set the final library name
SET( AMP_INC AMP )              # Set an optional subfolder for includes (e.g. include/name/...)
IF ( NOT ${PROJ}_INSTALL_DIR )
    SET( ${PROJ}_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}" )
ENDIF()


# Initialize the project (should disable langagues that the TPL builder will enable)
PROJECT( ${PROJ} LANGUAGES )


# Prevent users from building in place
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
    MESSAGE( FATAL_ERROR "Building code in place is a bad idea" )
ENDIF()


# Set source/install paths
SET( ${PROJ}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )
SET( ${PROJ}_BUILD_DIR  "${CMAKE_CURRENT_BINARY_DIR}" )
INCLUDE_DIRECTORIES( "${${PROJ}_INSTALL_DIR}/include" )


# Disable link-time optimization
SET( DISABLE_LTO TRUE)


# Build the TPLs if necessary
IF ( TPL_URL )
    INCLUDE( cmake/configure_TPLs.cmake )
    CONFIGURE_TPL_BUILDER()
ENDIF()


# Load the TPL package
IF ( NOT ONLY_BUILD_DOCS )
    IF ( NOT TPL_DIRECTORY )
        MESSAGE( FATAL_ERROR "TPL_DIRECTORY must be specified")
    ENDIF()
    SET( CMAKE_MODULE_PATH "${TPL_DIRECTORY}" ${CMAKE_MODULE_PATH} )
    IF ( NOT DEFINED TPL_LIST_REQUIRED )
        SET( TPL_LIST_REQUIRED LAPACK LAPACK_WRAPPERS STACKTRACE )
    ENDIF()
    IF ( NOT DEFINED TPL_LIST_OPTIONAL )
        SET( TPL_LIST_OPTIONAL MPI BOOST ZLIB EIGEN PETSC HDF5 SILO HYPRE TRILINOS SUNDIALS LIBMESH X11 OGRE TIMER SAMRAI )
    ENDIF()
    STRING( REGEX REPLACE  ","  ";"  TPL_LIST_REQUIRED  "${TPL_LIST_REQUIRED}" )
    STRING( REGEX REPLACE  ","  ";"  TPL_LIST_OPTIONAL  "${TPL_LIST_OPTIONAL}" )
    FIND_PACKAGE( TPLs REQUIRED ${TPL_LIST_REQUIRED} OPTIONAL_COMPONENTS ${TPL_LIST_OPTIONAL} )
    MESSAGE("TPL Found")
    SET( TPL_LIST_FOUND )
    FOREACH( tmp ${TPL_LIST} )
        IF ( TPL_FOUND_${tmp} )
            SET( TPL_LIST_FOUND ${TPL_LIST_FOUND} ${tmp} )
            SET( USE_EXT_${tmp} TRUE )
            ADD_DEFINITIONS( -DUSE_EXT_${tmp} )
        ENDIF()
    ENDFOREACH()
    MESSAGE("  TPLs Required: ${TPL_LIST_REQUIRED}")
    MESSAGE("  TPLs Supported: ${TPL_LIST_OPTIONAL}")
    MESSAGE("  TPLs Built: ${TPL_LIST}")
    MESSAGE("  TPLs Included: ${TPL_LIST_FOUND}")
    MESSAGE("  Include Paths: ${TPL_INCLUDE_DIRS}")
    MESSAGE("  Libraries: ${TPL_LIBRARIES}")
ENDIF()


# Include the TPL inclde paths and libraries
INCLUDE_DIRECTORIES( ${TPL_INCLUDE_DIRS} )
SET( TPL_LIBS ${TPL_LIBRARIES} )


# Set testing paramaters
SET( DROP_METHOD "http" )
SET( DROP_SITE "" )
SET( DROP_LOCATION "/CDash/submit.php?project=AMP" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
ENABLE_TESTING()
INCLUDE( CTest )


# Check for CUDA
CHECK_ENABLE_FLAG( USE_CUDA 0 )
IF ( USE_CUDA )
    ADD_DEFINITIONS( -DUSE_CUDA )
    ENABLE_LANGUAGE( CUDA )
    SET( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CXX_STD_FLAG}" )
ENDIF()


# Configure external libraries that are not configured through the TPL builder
IF ( NOT ONLY_BUILD_DOCS )
    INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake" )
    CONFIGURE_LINE_COVERAGE()
    CONFIGURE_OPENMP()
    IF ( NOT TPL_FOUND_TIMER )
        INCLUDE( "${TPL_DIRECTORY}/cmake/FindTimer.cmake" )
        CONFIGURE_TIMER( FALSE "${AMP_INSTALL_DIR}/include" FALSE )
    ENDIF()
    INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SharedPtr.cmake" )
    CONFIGURE_SHARED_PTR( "${${PROJ}_INSTALL_DIR}/include/AMP/utils" AMP )
    # Configure AMP (all external libraries must be configured first)
    CONFIGURE_AMP()
ENDIF()


# Create the target for documentation
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
CHECK_ENABLE_FLAG( USE_EXT_DOXYGEN 1 )
CHECK_ENABLE_FLAG( USE_LATEX 1 )
FILE( MAKE_DIRECTORY "${${PROJ}_INSTALL_DIR}/doc" )
IF ( USE_EXT_DOXYGEN )
    SET( DOXYFILE_IN "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/Doxyfile.in" )
    SET( DOXY_HEADER_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/html/header.html" )
    SET( DOXY_FOOTER_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/html/footer.html" )
    SET( DOXY_LATEX_HEADER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/doxygen/header.tex )
    SET( DOXYFILE_OUTPUT_DIR "${${PROJ}_INSTALL_DIR}/doc" )
    SET( DOXYFILE_SRC_HTML_DIR "${CMAKE_CURRENT_SOURCE_DIR}/doxygen/html" )
    SET( DOXYFILE_SOURCE_DIR "${${PROJ}_SOURCE_DIR}" )
    SET( DOXYFILE_EXTRA_SOURCES CACHE INTERNAL "doxyfile_extra_sources" FORCE )
    SET( REL_PACKAGE_HTML "" )
    SET( DOXYGEN_MACROS "USE_AMP_UTILS USE_AMP_MESH USE_AMP_DISCRETIZATION USE_AMP_VECTORS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_MATRICES USE_AMP_MATERIALS USE_AMP_OPERATORS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_SOLVERS USE_AMP_TIME_INTEGRATORS USE_AMP_GRAPHICS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_MPI USE_EXT_LAPACK USE_EXT_BLAS USE_EXT_PETSC" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_LIBMESH USE_EXT_SUNDIALS USE_EXT_TRILINOS" )
    SET(DOXYFILE_LATEX "NO")
    INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/UseDoxygen.cmake )
    IF ( NOT DOXYGEN_FOUND )
        MESSAGE( WARNING "Doxygen not found, disabling documentation.  Set USE_EXT_DOXYGEN=0 to disable this message")
        SET( USE_EXT_DOXYGEN OFF )
    ELSE()
        MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
    ENDIF()
ENDIF()


# Create custom targets for build-test, check, and distclean
SET( EXCLUDE_TESTS_FROM_ALL 0 )
ADD_CUSTOM_TARGET( build-test )
ADD_CUSTOM_TARGET( build-examples )
ADD_CUSTOM_TARGET( check COMMAND  make test  )
ADD_DISTCLEAN( src ampdir AMP compile_commands.json AMP-Data.tar.gz tpl-build tpl-builder )


# Add a timer for the build
CHECK_ENABLE_FLAG( TIME_BUILD 0 )
IF ( TIME_BUILD )
    SET_PROPERTY(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
    MESSAGE( STATUS "Adding build timers" )
ENDIF()


# Add the src directory
IF ( NOT ONLY_BUILD_DOCS )
    WRITE_REPO_VERSION()
    SET( EXCLUDE_TESTS_FROM_ALL 1 )
    ADD_SUBDIRECTORY( src )
ENDIF()


# Add the cppcheck tests
SET( CPPCHECK_TIMEOUT 700 )
SET( CPPCHECK_INCLUDE "-I${TIMER_INCLUDE}" "-I${${PROJ}_INSTALL_DIR}/include" )
FIND_PACKAGE( Cppcheck )


# Add the cppclean tests
SET( CPPCLEAN_OPTIONS )
SET( CPPCLEAN_EXCLUDE 
    libmeshGenerators.h meshGenerators.h meshTests.h
    DOFManager_tests.h test_Discretization.h
    Vector.h VectorDataIterator.h PetscVector.h NativePetscVector.inline.h
    Material.h
)
SET( CPPCLEAN_SUPPRESSIONS 
    ".inline.h' does not need to be #included"
    "tmpl.h' does not need to be #included"
    ".I' does not need to be #included"
    "'ProfilerApp.h' does not need to be #included"
    "'utils/shared_ptr.h' does not need to be #included"
    "'Utilities.h' does not need to be #included"
    "'UtilityMacros.h' does not need to be #included"
    "'utils/Counter.h' does not need to be #included"
    "'libmeshGenerators.h' does not need to be #included"
    "should #include header file '.*/src/.*'"
)
SET( CPPCLEAN_UNNECESSARY_INCLUDE 1 )
SET( CPPCLEAN_EXTRA_INCLUDE 1 )
SET( CPPCLEAN_SHOULD_INCLUDE 1 )
SET( CPPCLEAN_INCLUDE_NOT_FOUND 1 )
SET( CPPCLEAN_FUN_NOT_FOUND 0 )
SET( CPPCLEAN_DECLARED 0 )
SET( CPPCLEAN_STATIC 0 )
SET( CPPCLEAN_UNUSED_VARIABLE 0 )
SET( CPPCLEAN_UNKNOWN 1 )
SET( CPPCLEAN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}" )
FIND_PACKAGE( Cppclean )


# Create release target
IF ( AMP_MAJOR_VERSION AND AMP_MINOR_VERSION )
    ADD_CUSTOM_TARGET( release )
    INCLUDE( "${${PROJ}_INSTALL_DIR}/${PROJ}_Version.cmake" )
    SET( AMP_NAME AMP-${AMP_MAJOR_VERSION}.${AMP_MINOR_VERSION}.${AMP_BUILD_VERSION} )
    SET( RELEASE_FILENAME ${AMP_NAME}.tar )
    SET( TMP_DIR "${CMAKE_CURRENT_BINARY_DIR}/tmp" )
    FILE(WRITE  "${TMP_DIR}/release.cmake" "# Create the release tar\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory \"${CMAKE_CURRENT_SOURCE_DIR}\" ${AMP_NAME})\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${AMP_NAME}/.hg )\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E remove -f ${AMP_NAME}/.hgtags )\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E remove -f ${AMP_NAME}/.hgignore )\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E copy version.cmake ${AMP_NAME}/AMP_Version.cmake)\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E tar -cf ${AMP_NAME}.tar ${AMP_NAME})\n" )
    FILE(APPEND "${TMP_DIR}/release.cmake" "execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${AMP_NAME}.tar \"${AMP_INSTALL_DIR}/${${AMP_NAME}}\" )\n" )
    ADD_CUSTOM_COMMAND( TARGET release PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -P "${TMP_DIR}/release.cmake"
        WORKING_DIRECTORY "${TMP_DIR}" 
    )
ENDIF()


# Install the documentation
IF ( USE_EXT_DOXYGEN )
    ADD_DEPENDENCIES( doxygen latex_docs )
    ADD_DEPENDENCIES( doc latex_docs doxygen )
ENDIF()


# Save the configuration info for applications 
IF ( NOT ONLY_BUILD_DOCS )
    INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/save_variables.cmake )
    CONFIGURE_FILE( cmake/FindAMP.template.cmake "${${PROJ}_INSTALL_DIR}/FindAMP.cmake" @ONLY )
    COPY_DATA_FILE( "${TPL_MACRO_CMAKE}" "${${PROJ}_INSTALL_DIR}/cmake/macros.cmake" )
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SharedPtr.cmake" "${${PROJ}_INSTALL_DIR}/cmake/SharedPtr.cmake" )
    SAVE_CMAKE_FLAGS( "${${PROJ}_INSTALL_DIR}/amp.cmake" )
ENDIF()


