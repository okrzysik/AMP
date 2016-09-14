# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
CMAKE_POLICY(SET CMP0004 OLD)
IF ( ${CMAKE_MAJOR_VERSION} EQUAL 3 )
    CMAKE_POLICY(SET CMP0026 OLD)
    CMAKE_POLICY(SET CMP0053 OLD)
ENDIF()


MESSAGE ("===============")
MESSAGE ("Configuring AMP")
MESSAGE ("===============")


# Set the project name
SET( PROJ AMP )
SET( AMP_LIB amp )


# Initialize the project (should disable langagues that the TPL builder will enable)
PROJECT( ${PROJ} LANGUAGES )


# Prevent users from building in place
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
    MESSAGE( FATAL_ERROR "Building code in place is a bad idea" )
ENDIF()


# Set source/install paths
SET( ${PROJ}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src" )
SET( ${PROJ}_BUILD_DIR  "${CMAKE_CURRENT_BINARY_DIR}" )
IF( ${PROJ}_INSTALL_DIR )
    SET( ${PROJ}_INSTALL_DIR "${${PROJ}_INSTALL_DIR}" )
ELSEIF( NOT ${PROJ}_INSTALL_DIR )
    SET( ${PROJ}_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/ampdir" )
ENDIF()
INCLUDE_DIRECTORIES( "${${PROJ}_INSTALL_DIR}/include" )


# Load the TPL package
IF ( NOT ONLY_BUILD_DOCS )
    IF ( NOT TPL_DIRECTORY )
        MESSAGE( FATAL_ERROR "TPL_DIRECTORY must be specified")
    ENDIF()
    SET( CMAKE_MODULE_PATH "${TPL_DIRECTORY}" ${CMAKE_MODULE_PATH} )
    IF ( NOT DEFINED TPL_LIST_REQUIRED )
        SET( TPL_LIST_REQUIRED LAPACK )
    ENDIF()
    IF ( NOT DEFINED TPL_LIST_OPTIONAL )
        SET( TPL_LIST_OPTIONAL MPI BOOST ZLIB EIGEN PETSC HDF5 HYPRE TRILINOS SUNDIALS LIBMESH X11 OGRE TIMER )
    ENDIF()
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


# Configure external libraries that are not configured through the TPL builder
IF ( NOT ONLY_BUILD_DOCS )
    INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake" )
    CONFIGURE_LINE_COVERAGE()
    CONFIGURE_CUDA()
    IF ( NOT TPL_FOUND_TIMER )
        INCLUDE( "${TPL_DIRECTORY}/cmake/FindTimer.cmake" )
        CONFIGURE_TIMER( FALSE "${AMP_INSTALL_DIR}/include" )
    ENDIF()
    INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SharedPtr.cmake" )
    CONFIGURE_SHARED_PTR( "${${PROJ}_INSTALL_DIR}/include/utils" AMP )
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
ADD_DISTCLEAN( src ampdir AMP compile_commands.json )
MACRO( ADD_AMP_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_AMP_TEST ( ${EXENAME} ${ARGN} )
    ADD_AMP_TEST_PARALLEL ( ${EXENAME} 2 ${ARGN} )
    ADD_AMP_TEST_PARALLEL ( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Add the src directory
IF ( NOT ONLY_BUILD_DOCS )
    SET( EXCLUDE_TESTS_FROM_ALL 1 )
    ADD_SUBDIRECTORY( src )
ENDIF()


# Add the cppcheck tests
SET( CPPCHECK_INCLUDE "-I${TIMER_INCLUDE}" )
SET( CPPCHECK_OPTIONS -q --enable=all --force --suppress=missingIncludeSystem 
    "--suppressions-list=${CMAKE_CURRENT_SOURCE_DIR}/cppcheckSuppressionFile" )
IF ( CXX_STD STREQUAL 98 )
    SET( CPPCHECK_OPTIONS ${CPPCHECK_OPTIONS} --std=c99 --std=c++03 --std=posix )
ELSEIF ( CXX_STD STREQUAL 11 )
    SET( CPPCHECK_OPTIONS ${CPPCHECK_OPTIONS} --std=c11 --std=c++11 --std=posix )
ELSEIF ( CXX_STD STREQUAL 14 )
    SET( CPPCHECK_OPTIONS ${CPPCHECK_OPTIONS} --std=c11 --std=c++11 --std=posix )
ENDIF()
SET( CPPCHECK_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}" )
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
    "'utils/Utilities.h' does not need to be #included"
    "'utils/Counter.h' does not need to be #included"
)
SET( CPPCLEAN_UNNECESSARY_INCLUDE 1 )
SET( CPPCLEAN_EXTRA_INCLUDE 1 )
SET( CPPCLEAN_SHOULD_INCLUDE 0 )
SET( CPPCLEAN_INCLUDE_NOT_FOUND 0 )
SET( CPPCLEAN_FUN_NOT_FOUND 0 )
SET( CPPCLEAN_DECLARED 0 )
SET( CPPCLEAN_STATIC 0 )
SET( CPPCLEAN_UNUSED_VARIABLE 0 )
SET( CPPCLEAN_UNKNOWN 1 )
SET( CPPCLEAN_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}" )
FIND_PACKAGE( Cppclean )


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
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteRepoVersion.cmake" "${${PROJ}_INSTALL_DIR}/cmake/WriteRepoVersion.cmake" )
    SAVE_CMAKE_FLAGS( "${${PROJ}_INSTALL_DIR}/amp.cmake" )
ENDIF()


