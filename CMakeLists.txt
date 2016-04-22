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


# Check if we used the TPL builder (required before PROJECT)
SET( PROJ AMP )
SET( AMP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
SET( AMP_BUILD_DIR  "${CMAKE_CURRENT_BINARY_DIR}" )
IF( PREFIX )
    SET( AMP_INSTALL_DIR "${PREFIX}" )
ELSEIF( NOT AMP_INSTALL_DIR )
    SET( AMP_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/ampdir" )
ENDIF()
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake" )
IF ( TPL_DIRECTORY )
    INCLUDE( "${TPL_DIRECTORY}/TPLs.cmake" )
ENDIF()


# Set the project name
PROJECT (AMP)


# Check if we are only compiling docs
CHECK_ENABLE_FLAG( ONLY_BUILD_DOCS 0 )


# Enable the required compilers
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)


# Re-include TPLs.cmake to fix bug with CMAKE_COMPILE_MODE being reset
IF ( TPL_DIRECTORY )
    INCLUDE( ${TPL_DIRECTORY}/TPLs.cmake )
    CONFIGURE_TPLs()
ENDIF()


# Call TPLsList to get the TPLs
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/TPLsList.cmake" )
MESSAGE("TPL Location : ${TPL_DIRECTORY}")
MESSAGE("COMPILE_MODE : ${CMAKE_BUILD_TYPE} " )


# Call PackageList to get the packages, configure the directories and compilers, etc
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/PackagesList.cmake" )


# Set testing paramaters
SET( DROP_METHOD "http" )
SET( DROP_SITE "" )
SET( DROP_LOCATION "/CDash/submit.php?project=AMP" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
ENABLE_TESTING()
INCLUDE( CTest )


# Configure external libraries
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_MPI()     # MPI must be before other libraries
    CONFIGURE_EIGEN_LIBRARIES()
    CONFIGURE_ARMADILLO_LIBRARIES()
    CONFIGURE_TRILINOS_LIBRARIES()
    CONFIGURE_SUNDIALS_LIBRARIES()
    CONFIGURE_PETSC_LIBRARIES()
    CONFIGURE_HYPRE_LIBRARIES()
    CONFIGURE_X11_LIBRARIES()
    CONFIGURE_BLAS_AND_LAPACK()
    CONFIGURE_LIBMESH()
    CONFIGURE_NEK()     # should go before MOAB since it requires MOAB.
    CONFIGURE_MOAB()
    CONFIGURE_DENDRO()
    CONFIGURE_SILO()
    CONFIGURE_HDF5()
    CONFIGURE_NETCDF()
    CONFIGURE_BOOST()
    CONFIGURE_ZLIB()
    CHECK_ENABLE_FLAG( USE_EXT_TIMER 1 )
    CONFIGURE_TIMER( ${USE_EXT_TIMER} "${AMP_INSTALL_DIR}/include" )
    CONFIGURE_SHARED_PTR( "${AMP_INSTALL_DIR}/include/utils" AMP )
ENDIF()


# Configure AMP (all external libraries must be configured first)
CONFIGURE_AMP()


# Create the target for documentation
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
CHECK_ENABLE_FLAG( USE_EXT_DOXYGEN 1 )
CHECK_ENABLE_FLAG( USE_LATEX 1 )
FILE( MAKE_DIRECTORY "${AMP_INSTALL_DIR}/doc" )
IF ( USE_EXT_DOXYGEN )
    SET( DOXYFILE_IN ${AMP_SOURCE_DIR}/doxygen/Doxyfile.in)
    SET( DOXY_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/html/header.html )
    SET( DOXY_FOOTER_FILE ${AMP_SOURCE_DIR}/doxygen/html/footer.html )
    SET( DOXY_LATEX_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/header.tex )
    SET( DOXYFILE_OUTPUT_DIR "${AMP_INSTALL_DIR}/doc" )
    SET( DOXYFILE_SRC_HTML_DIR "${AMP_SOURCE_DIR}/doxygen/html" )
    SET( DOXYFILE_SOURCE_DIR "${AMP_SOURCE_DIR}/src" )
    SET( DOXYFILE_EXTRA_SOURCES CACHE INTERNAL "doxyfile_extra_sources" FORCE )
    SET( REL_PACKAGE_HTML "" )
    SET( DOXYGEN_MACROS "USE_AMP_UTILS USE_AMP_MESH USE_AMP_DISCRETIZATION USE_AMP_VECTORS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_MATRICES USE_AMP_MATERIALS USE_AMP_OPERATORS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_SOLVERS USE_AMP_TIME_INTEGRATORS" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_MPI USE_EXT_LAPACK USE_EXT_BLAS USE_EXT_PETSC" )
    SET( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_LIBMESH USE_EXT_SUNDIALS USE_EXT_TRILINOS" )
    SET(DOXYFILE_LATEX "NO")
    INCLUDE( ${AMP_SOURCE_DIR}/cmake/UseDoxygen.cmake )
    IF ( NOT DOXYGEN_FOUND )
        MESSAGE( WARNING "Doxygen not found, disabling documentation.  Set USE_EXT_DOXYGEN=0 to disable this message")
        SET( USE_EXT_DOXYGEN OFF )
    ELSE()
        MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
    ENDIF()
ENDIF()


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
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppcheck.cmake" )


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
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppclean.cmake" )


# Install the documentation
IF ( USE_EXT_DOXYGEN )
    ADD_DEPENDENCIES( doxygen latex_docs )
    ADD_DEPENDENCIES( doc latex_docs doxygen )
ENDIF()


# Save the configuration info for applications 
IF ( NOT ONLY_BUILD_DOCS )
    INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/save_variables.cmake )
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake" "${${PROJ}_INSTALL_DIR}/cmake/macros.cmake" )
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cppcheck.cmake" "${${PROJ}_INSTALL_DIR}/cmake/cppcheck.cmake" )
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/SharedPtr.cmake" "${${PROJ}_INSTALL_DIR}/cmake/SharedPtr.cmake" )
    COPY_DATA_FILE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/WriteRepoVersion.cmake" "${${PROJ}_INSTALL_DIR}/cmake/WriteRepoVersion.cmake" )
    SAVE_CMAKE_FLAGS( "${${PROJ}_INSTALL_DIR}/amp.cmake" )
ENDIF()


