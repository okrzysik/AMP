# Set some CMake properties
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
CMAKE_POLICY(SET CMP0014 OLD)


MESSAGE ("===============")
MESSAGE ("Configuring AMP")
MESSAGE ("===============")


# Prevent users from building in place
IF ("${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}" )
  MESSAGE ( FATAL_ERROR "Building AMP in place is a bad idea" )
ENDIF ()


PROJECT (AMP)
ENABLE_LANGUAGE (Fortran)

# Set testing paramaters
SET ( DROP_METHOD "http" )
SET ( DROP_SITE "thyme.ornl.gov" )
SET ( DROP_LOCATION "/CDash/submit.php?project=AMP" )
SET ( TRIGGER_SITE "" )
SET ( DROP_SITE_CDASH TRUE )
ENABLE_TESTING ()
INCLUDE (CTest)

# Identify the source and build directories
SET (AMP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
SET (AMP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
IF ( NOT AMP_INSTALL_DIR )
  SET ( AMP_INSTALL_DIR ${AMP_BUILD_DIR}/ampdir )
ENDIF()
MESSAGE ( "Installing AMP in "${AMP_INSTALL_DIR} )
SET (CMAKE_MODULE_PATH ${AMP_SOURCE_DIR} ${AMP_SOURCE_DIR}/cmake)

# Set Cmake include files
INCLUDE ( CheckIncludeFile )
INCLUDE ( macros )
INCLUDE ( libraries )


# Create custom targets for build-test, check, and distclean
ADD_CUSTOM_TARGET ( build-test )
ADD_CUSTOM_TARGET ( check COMMAND  make test  )
ADD_DISTCLEAN ()


# Check the compile mode and compile flags
VERIFY_VARIABLE ( "COMPILE_MODE" )
IF ( ${COMPILE_MODE} STREQUAL "debug" )
  SET_DEBUG_MACROS()
ELSE ()
  IF ( ${COMPILE_MODE} STREQUAL "optimized" )
    SET_OPTIMIZED_MACROS()
  ELSE ()
    MESSAGE ( FATAL_ERROR "COMPILE_MODE must be either debug or optimized" )
  ENDIF()
ENDIF()


# Configure external libraries
CONFIGURE_MPI ()   # MPI must be before other libraries
CONFIGURE_TRILINOS_LIBRARIES ()
CONFIGURE_SUNDIALS_LIBRARIES ()
CONFIGURE_PETSC_LIBRARIES ()
CONFIGURE_HYPRE_LIBRARIES ()
CONFIGURE_X11_LIBRARIES ()
CONFIGURE_BLAS ()
CONFIGURE_LAPACK ()
CONFIGURE_LIBMESH ()
CONFIGURE_NEK ()  # should go before MOAB since it requires MOAB.
CONFIGURE_MOAB ()
CONFIGURE_SILO ()
CONFIGURE_HDF5 ()
CONFIGURE_BOOST ()


# Configure AMP (all external libraries must be configured first)
CONFIGURE_AMP ()


# Some more configure options
CONFIGURE_TIMERS ()
CONFIGURE_LINE_COVERAGE ()
ADD_DEFINITIONS ( -DCMAKE_CONFIGURED )
INCLUDE (FortranCInterface)
FORTRANCINTERFACE_HEADER ( ${AMP_INSTALL_DIR}/include/utils/FC.h MACRO_NAMESPACE "FC_" )
if( (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") AND (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") )
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-common")
	set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-common")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-common")
endif( (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") AND (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") )


# Create the target for documentation
SET ( DOXYFILE_IN ${AMP_SOURCE_DIR}/doxygen/Doxyfile.in)
SET ( DOXY_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/html/header.html )
SET ( DOXY_FOOTER_FILE ${AMP_SOURCE_DIR}/doxygen/html/footer.html )
SET ( DOXY_LATEX_HEADER_FILE ${AMP_SOURCE_DIR}/doxygen/header.tex )
SET ( DOXYFILE_OUTPUT_DIR ${AMP_INSTALL_DIR}/doc/ )
SET ( DOXYFILE_SRC_HTML_DIR ${AMP_SOURCE_DIR}/doxygen/html )
SET ( DOXYFILE_SOURCE_DIR ${AMP_SOURCE_DIR}/src )
SET ( REL_PACKAGE_HTML "" )
SET ( DOXYGEN_MACROS "USE_SUNDIALS:=1" )
INCLUDE (UseDoxygen)


# Add the src directory
ADD_SUBDIRECTORY ( src )


# Save the configuration info for applications 
SAVE_CMAKE_FLAGS()


