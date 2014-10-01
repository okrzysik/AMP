# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#CMAKE_POLICY(SET CMP0011 OLD)
#CMAKE_POLICY(SET CMP0014 OLD)
CMAKE_POLICY(SET CMP0004 OLD)

MESSAGE ("===============")
MESSAGE ("Configuring AMP")
MESSAGE ("===============")


# Check if we used the TPL builder (required before PROJECT)
SET( AMP_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libraries.cmake" )
IF ( TPL_DIRECTORY )
    INCLUDE( ${TPL_DIRECTORY}/TPLs.cmake )
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
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/TPLsList.cmake )


MESSAGE("TPL Location : ${TPL_DIRECTORY}")
MESSAGE("COMPILE_MODE : ${CMAKE_BUILD_TYPE} " )
# Call PackageList to get the pacakages, configure the directories and compilers, etc
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/PackagesList.cmake )


# Set testing paramaters
SET( DROP_METHOD "http" )
SET( DROP_SITE "thyme.ornl.gov" )
SET( DROP_LOCATION "/CDash/submit.php?project=AMP" )
SET( TRIGGER_SITE "" )
SET( DROP_SITE_CDASH TRUE )
ENABLE_TESTING()
INCLUDE(CTest)


# Configure external libraries
IF ( NOT ONLY_BUILD_DOCS )
    CONFIGURE_MPI()     # MPI must be before other libraries
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
    CONFIGURE_ZLIB()
    CONFIGURE_BOOST()
ENDIF()


# Configure AMP (all external libraries must be configured first)
CONFIGURE_AMP()


# Create the target for documentation
ADD_CUSTOM_TARGET( doc )
ADD_CUSTOM_TARGET( latex_docs )
CHECK_ENABLE_FLAG(USE_EXT_DOXYGEN 1 )
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
    MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
ENDIF()


# Add the src directory
IF ( NOT ONLY_BUILD_DOCS )
    SET( EXCLUDE_TESTS_FROM_ALL 1 )
    ADD_SUBDIRECTORY( src )
ENDIF()


# Install the documentation
IF ( USE_EXT_DOXYGEN )
    INCLUDE( ${AMP_SOURCE_DIR}/cmake/UseDoxygen.cmake )
    SET(DOXYFILE_LATEX "YES")
    ADD_DEPENDENCIES( doxygen latex_docs )
    ADD_DEPENDENCIES( doc latex_docs doxygen )
ENDIF()


# Save the configuration info for applications 
IF ( NOT ONLY_BUILD_DOCS )
    SAVE_CMAKE_FLAGS()
ENDIF()


