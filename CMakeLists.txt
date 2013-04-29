# Set some CMake properties    
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
CMAKE_POLICY(SET CMP0011 OLD)
CMAKE_POLICY(SET CMP0014 OLD)


MESSAGE ("===============")
MESSAGE ("Configuring AMP")
MESSAGE ("===============")


ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)


# Set the project name
PROJECT (AMP)


# Call TPLsList to get the TPLs
INCLUDE( ${CMAKE_CURRENT_SOURCE_DIR}/TPLsList.cmake )


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


# Set Cmake include files
INCLUDE( ${AMP_SOURCE_DIR}/cmake/macros.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/libraries.cmake )


# Configure external libraries
CONFIGURE_MPI()     # MPI must be before other libraries
CONFIGURE_TRILINOS_LIBRARIES()
CONFIGURE_SUNDIALS_LIBRARIES()
CONFIGURE_PETSC_LIBRARIES()
CONFIGURE_HYPRE_LIBRARIES()
CONFIGURE_X11_LIBRARIES()
CONFIGURE_BLAS()
CONFIGURE_LAPACK()
CONFIGURE_LIBMESH()
CONFIGURE_NEK()     # should go before MOAB since it requires MOAB.
CONFIGURE_MOAB()
CONFIGURE_DENDRO()
CONFIGURE_SILO()
CONFIGURE_HDF5()
CONFIGURE_NETCDF()
CONFIGURE_BOOST()


# Configure AMP (all external libraries must be configured first)
CONFIGURE_AMP()


# Create the target for documentation
CHECK_ENABLE_FLAG(USE_EXT_DOXYGEN 1 )
IF ( USE_EXT_DOXYGEN )
    SET ( DOXYFILE_IN ${AMP_SOURCE_DIR}/doxygen/Doxyfile.in)
    SET ( DOXY_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/html/header.html )
    SET ( DOXY_FOOTER_FILE ${AMP_SOURCE_DIR}/doxygen/html/footer.html )
    SET ( DOXY_LATEX_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/header.tex )
    SET ( DOXYFILE_OUTPUT_DIR ${AMP_INSTALL_DIR}/doc/ )
    SET ( DOXYFILE_SRC_HTML_DIR ${AMP_SOURCE_DIR}/doxygen/html )
    SET ( DOXYFILE_SOURCE_DIR ${AMP_SOURCE_DIR}/src )
    SET ( REL_PACKAGE_HTML "" )
    SET ( DOXYGEN_MACROS "USE_AMP_UTILS USE_AMP_MESH USE_AMP_DISCRETIZATION USE_AMP_VECTORS" )
    SET ( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_MATRICES USE_AMP_MATERIALS USE_AMP_OPERATORS" )
    SET ( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_AMP_SOLVERS USE_AMP_TIME_INTEGRATORS" )
    SET ( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_MPI USE_EXT_LAPACK USE_EXT_BLAS USE_EXT_PETSC" )
    SET ( DOXYGEN_MACROS "${DOXYGEN_MACROS} USE_EXT_LIBMESH USE_EXT_SUNDIALS USE_EXT_TRILINOS" )
    MESSAGE("DOXYGEN_MACROS = ${DOXYGEN_MACROS}")
    INCLUDE( ${AMP_SOURCE_DIR}/cmake/UseDoxygen.cmake )
ENDIF()


# Add the src directory
SET( EXCLUDE_TESTS_FROM_ALL 1 )
ADD_SUBDIRECTORY ( src )


# Save the configuration info for applications 
SAVE_CMAKE_FLAGS()


