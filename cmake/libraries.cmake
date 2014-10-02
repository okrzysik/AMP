INCLUDE( ${AMP_SOURCE_DIR}/cmake/Find_BLAS_LAPACK.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/FindPetsc.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/FindTrilinos.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/FindLibmesh.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/FindSundials.cmake )
INCLUDE( ${AMP_SOURCE_DIR}/cmake/configureAMP.cmake )
INCLUDE( CheckIncludeFile )


MACRO ( CONFIGURE_LINE_COVERAGE )
    SET ( COVERAGE_LIBS )
    IF ( USE_EXT_GCOV )
        SET ( COVERAGE_LIBS ${COVERAGE_LIBS} -lgcov )
    ENDIF ()
    IF ( ENABLE_GCOV )
        ADD_DEFINITIONS ( -fprofile-arcs -ftest-coverage )
        SET ( COVERAGE_LIBS ${COVERAGE_LIBS} -fprofile-arcs )
    ENDIF ()
ENDMACRO ()


MACRO ( CONFIGURE_TIMERS )
  CHECK_INCLUDE_FILE ( sys/times.h HAVE_SYS_TIMES_H )
  CHECK_INCLUDE_FILE ( windows.h HAVE_WINDOWS_H )
ENDMACRO ()


# Macro to configure TPLs compiled through the TPL builder
MACRO( CONFIGURE_TPLs )
    CHECK_ENABLE_FLAG(USE_EXT_BOOST 0 )
    IF ( BOOST_INSTALL_DIR )
        SET( USE_EXT_BOOST true )
        SET( BOOST_DIRECTORY "${BOOST_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_ZLIB 0 )
    IF ( ZLIB_INSTALL_DIR )
        SET( USE_EXT_ZLIB true )
        SET( ZLIB_DIRECTORY "${ZLIB_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_PETSC 0 )
    IF ( PETSC_INSTALL_DIR )
        SET( USE_EXT_PETSC true )
        SET( PETSC_DIRECTORY "${PETSC_INSTALL_DIR}" )
        SET( PETSC_ARCH "" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_HDF5 0 )
    IF ( HDF5_INSTALL_DIR )
        SET( USE_EXT_HDF5 true )
        SET( HDF5_DIRECTORY "${HDF5_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_SILO 0 )
    IF ( SILO_INSTALL_DIR )
        SET( USE_EXT_SILO true )
        SET( SILO_DIRECTORY "${SILO_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_HYPRE 0 )
    IF ( HYPRE_INSTALL_DIR )
        SET( USE_EXT_HYPRE true )
        SET( HYPRE_DIRECTORY "${HYPRE_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_LAPACK 0 )
    IF ( LAPACK_INSTALL_DIR )
        SET( USE_EXT_LAPACK true )
        SET( LAPACK_DIRECTORY "${LAPACK_INSTALL_DIR}" )
        SET( BLAS_LIBRARIES "${BLAS_LIBS}" )
        SET( LAPACK_LIBRARIES "${LAPACK_LIBS}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_TRILINOS 0 )
    IF ( TRILINOS_INSTALL_DIR )
        SET( USE_EXT_TRILINOS true )
        SET( TRILINOS_DIRECTORY "${TRILINOS_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_SUNDIALS 0 )
    IF ( SUNDIALS_INSTALL_DIR )
        SET( USE_EXT_SUNDIALS true )
        SET( SUNDIALS_DIRECTORY "${SUNDIALS_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_LIBMESH 0 )
    IF ( LIBMESH_INSTALL_DIR )
        SET( USE_EXT_LIBMESH true )
        SET( LIBMESH_DIRECTORY "${LIBMESH_INSTALL_DIR}" )
    ENDIF()
    CHECK_ENABLE_FLAG(USE_EXT_X11 0 )
    IF ( X11_INSTALL_DIR )
        SET( USE_EXT_X11 true )
        SET( X11_DIRECTORY "${X11_INSTALL_DIR}" )
    ENDIF()
ENDMACRO()


# Macro to find and configure boost (we only need the headers)
MACRO ( CONFIGURE_BOOST )
    # Determine if we want to use boost
    CHECK_ENABLE_FLAG(USE_EXT_BOOST 1 )
    IF ( USE_EXT_BOOST )
        # Check if we specified the boost directory
        IF ( BOOST_DIRECTORY )
            VERIFY_PATH ( ${BOOST_DIRECTORY} )
            VERIFY_PATH ( ${BOOST_DIRECTORY}/include )
            SET ( BOOST_INCLUDE ${BOOST_DIRECTORY}/include )
        ELSE()
            # Check the default path for boost
            VERIFY_PATH ( ${AMP_SOURCE_DIR}/../external/boost/include )
            SET ( BOOST_INCLUDE ${AMP_SOURCE_DIR}/../external/boost/include )
        ENDIF()
        INCLUDE_DIRECTORIES ( ${BOOST_INCLUDE} )
        ADD_DEFINITIONS ( "-D USE_EXT_BOOST" )
        MESSAGE( "Using boost" )
    ELSE()
        MESSAGE( FATAL_ERROR "boost headers are necessary for AMP" )
    ENDIF()
ENDMACRO()


# Macro to find and configure netcdf
MACRO ( CONFIGURE_NETCDF )
    CHECK_ENABLE_FLAG(USE_EXT_NETCDF 0 )
    IF ( USE_EXT_NETCDF )
        IF ( NETCDF_DIRECTORY )
            VERIFY_PATH ( ${NETCDF_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${NETCDF_DIRECTORY}/include )
            SET ( NETCDF_INCLUDE ${NETCDF_DIRECTORY}/include )
            FIND_LIBRARY ( NETCDF_NETCDF_LIB    NAMES netcdf    PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( NETCDF_HDF5_LIB      NAMES hdf5      PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( NETCDF_HL_LIB        NAMES hl        PATHS ${NETCDF_DIRECTORY}/lib  NO_DEFAULT_PATH )
            IF ( (NOT NETCDF_NETCDF_LIB) OR (NOT NETCDF_HDF5_LIB) OR (NOT NETCDF_HL_LIB)  )
                MESSAGE( "${NETCDF_NETCDF_LIB}" )
                MESSAGE( "${NETCDF_HDF5_LIB}" )
                MESSAGE( "${NETCDF_HL_LIB}" )
            ENDIF()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for netcdf is not yet supported.  Use -D NETCDF_DIRECTORY=" )
        ENDIF()
        SET ( NETCDF_LIBS
            ${NETCDF_NETCDF_LIB} 
            ${NETCDF_HDF5_LIB} 
            ${NETCDF_HL_LIB} 
            ${NETCDF_NETCDF_LIB} 
            ${NETCDF_HDF5_LIB} 
            ${NETCDF_HL_LIB} 
        )
        ADD_DEFINITIONS ( "-D USE_NETCDF" )
        MESSAGE( "Using netcdf" )
        MESSAGE( "   ${NETCDF_LIBS}" )
    ENDIF()
ENDMACRO()


# Macro to find and configure the trilinos libraries
MACRO ( CONFIGURE_TRILINOS_LIBRARIES )
    # Determine if we want to use trilinos
    CHECK_ENABLE_FLAG(USE_EXT_TRILINOS 1 )
    IF ( USE_EXT_TRILINOS )
        # Check if we specified the trilinos directory
        IF ( TRILINOS_DIRECTORY )
            VERIFY_PATH( ${TRILINOS_DIRECTORY} )
            MESSAGE( " Trilinos Directory " ${TRILINOS_DIRECTORY} )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for trilinos is not yet supported.  Use -D TRILINOS_DIRECTORY=" )
        ENDIF()
        # Add the include directories
        INCLUDE_DIRECTORIES( ${TRILINOS_DIRECTORY}/include )
        SET ( TRILINOS_INCLUDE ${TRILINOS_DIRECTORY}/include )
        # Get the trilinos version
        TRILINOS_GET_VERSION()
        MESSAGE("Found Trilinos version ${TRILINOS_VERSION}")
        # Set the subpackages we want to use
        TRILINOS_SET_SUBPACKAGES()
        # Get the trilinos libraries
        IF (EXISTS "${TRILINOS_DIRECTORY}/lib/cmake/Trilinos/TrilinosConfig.cmake")
            INCLUDE("${TRILINOS_DIRECTORY}/lib/cmake/Trilinos/TrilinosConfig.cmake")
            SET( TRILINOS_LIBS )
            FOREACH (value2 ${Trilinos_LIBRARIES})
                UNSET( value CACHE )
                FIND_LIBRARY( value  NAMES ${value2}  PATHS ${Trilinos_LIBRARY_DIRS}  NO_DEFAULT_PATH )   
                SET( TRILINOS_LIBS ${TRILINOS_LIBS} ${value} )
            ENDFOREACH (value2)           
        ELSE()
            TRILINOS_SET_LIBRARIES()
        ENDIF()
        ADD_DEFINITIONS( "-D USE_EXT_TRILINOS" )
        MESSAGE( "Using trilinos" )
        MESSAGE( "   ${TRILINOS_LIBS}" )
    ELSE()
        MESSAGE("Configuring without trilinos")
    ENDIF()
ENDMACRO ()


# Macro to find and configure the silo libraries
MACRO ( CONFIGURE_SILO )
    # Determine if we want to use silo
    CHECK_ENABLE_FLAG(USE_EXT_SILO 1 )
    IF ( USE_EXT_SILO )
        # Check if we specified the silo directory
        IF ( SILO_DIRECTORY )
            VERIFY_PATH ( ${SILO_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${SILO_DIRECTORY}/include )
            SET ( SILO_INCLUDE ${SILO_DIRECTORY}/include )
            FIND_LIBRARY ( SILO_LIB  NAMES siloh5  PATHS ${SILO_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( "Default search for silo is not yet supported")
            MESSAGE( "Use -D SILO_DIRECTORY=" FATAL_ERROR)
        ENDIF()
        SET ( SILO_LIBS
            ${SILO_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_SILO" )  
        MESSAGE( "Using silo" )
        MESSAGE( "   ${SILO_LIB}" )
    ENDIF ()
ENDMACRO ()


# Macro to find and configure the hdf5 libraries
MACRO ( CONFIGURE_HDF5 )
    # Determine if we want to use hdf5
    CHECK_ENABLE_FLAG(USE_EXT_HDF5 1 )
    IF ( USE_EXT_HDF5 )
        # Check if we specified the silo directory
        IF ( HDF5_DIRECTORY )
            VERIFY_PATH ( ${HDF5_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${HDF5_DIRECTORY}/include )
            SET ( HDF5_INCLUDE ${HDF5_DIRECTORY}/include )
            FIND_LIBRARY ( HDF5_LIB    NAMES hdf5    PATHS ${HDF5_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HDF5_HL_LIB NAMES hdf5_hl PATHS ${HDF5_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for hdf5 is not yet supported.  Use -D HDF5_DIRECTORY=" )
        ENDIF()
        SET ( HDF5_LIBS
            ${HDF5_HL_LIB}
            ${HDF5_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_HDF5" )  
        MESSAGE( "Using hdf5" )
        MESSAGE( "   ${HDF5_LIB}" )
    ENDIF()
ENDMACRO ()

# Macro to find and configure zlib
MACRO ( CONFIGURE_ZLIB )
    # Determine if we want to use zlib
    CHECK_ENABLE_FLAG(USE_EXT_ZLIB 1 )
    IF ( USE_EXT_ZLIB )
        # Check if we specified the silo directory
        IF ( ZLIB_DIRECTORY )
            VERIFY_PATH ( ${ZLIB_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${ZLIB_DIRECTORY}/include )
            SET ( ZLIB_INCLUDE ${ZLIB_DIRECTORY}/include )
            FIND_LIBRARY ( ZLIB_LIB    NAMES z    PATHS ${ZLIB_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            FIND_LIBRARY ( ZLIB_LIB    NAMES z ) 
        ENDIF()
        SET ( ZLIB_LIBS
            ${ZLIB_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_ZLIB" )  
        MESSAGE( "Using zlib" )
        MESSAGE( "   ${ZLIB_LIB}" )
    ENDIF()
ENDMACRO ()

# Macro to find and configure the X11 libraries
MACRO ( CONFIGURE_X11_LIBRARIES )
    # Determine if we want to use X11
    CHECK_ENABLE_FLAG(USE_EXT_X11 1 )
    IF ( USE_EXT_X11 )
        # Check if we specified the silo directory
        IF ( X11_DIRECTORY )
            VERIFY_PATH ( ${X11_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${X11_DIRECTORY}/include )
            SET ( X11_INCLUDE ${X11_DIRECTORY}/include )
            FIND_LIBRARY ( X11_SM_LIB  NAMES SM  PATHS ${X11_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( X11_ICE_LIB NAMES ICE PATHS ${X11_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( X11_X11_LIB NAMES X11 PATHS ${X11_DIRECTORY}/lib  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for X11 is not yet supported.  Use -D X11_DIRECTORY=" )
        ENDIF()
        SET ( X11_LIBS
            ${X11_SM_LIB}
            ${X11_ICE_LIB}
            ${X11_X11_LIB} 
        )
        ADD_DEFINITIONS ( "-D USE_EXT_X11" )  
        MESSAGE( "Using X11" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure the MPI libraries
MACRO ( CONFIGURE_MPI )
    # Determine if we want to use MPI
    CHECK_ENABLE_FLAG(USE_EXT_MPI 1 )
    IF ( USE_EXT_MPI )
        # Check if we specified the MPI directory
        IF ( MPI_DIRECTORY )
            # Check the provided MPI directory for include files and the mpi executable
            VERIFY_PATH ( ${MPI_DIRECTORY} )
            SET ( MPI_INCLUDE_PATH ${MPI_DIRECTORY}/include )
            VERIFY_PATH ( ${MPI_INCLUDE_PATH} )
            IF ( NOT EXISTS ${MPI_INCLUDE_PATH}/mpi.h )
                MESSAGE( FATAL_ERROR "mpi.h not found in ${MPI_INCLUDE_PATH}/include" )
            ENDIF ()
            INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
            SET ( MPI_INCLUDE ${MPI_INCLUDE_PATH} )
            IF ( MPIEXEC ) 
                # User specified the MPI command directly, use as is
            ELSEIF ( MPIEXEC_CMD )
                # User specified the name of the MPI executable
                SET ( MPIEXEC ${MPI_DIRECTORY}/bin/${MPIEXEC_CMD} )
                IF ( NOT EXISTS ${MPIEXEC} )
                    MESSAGE( FATAL_ERROR "${MPIEXEC_CMD} not found in ${MPI_DIRECTORY}/bin" )
                ENDIF ()
            ELSE ()
                # Search for the MPI executable in the current directory
                FIND_PROGRAM ( MPIEXEC  NAMES mpiexec mpirun lamexec  PATHS ${MPI_DIRECTORY}/bin  NO_DEFAULT_PATH )
                IF ( NOT MPIEXEC )
                    MESSAGE( FATAL_ERROR "Could not locate mpi executable" )
                ENDIF()
            ENDIF ()
            # Set MPI flags
            IF ( NOT MPIEXEC_NUMPROC_FLAG )
                SET( MPIEXEC_NUMPROC_FLAG "-np" )
            ENDIF()
        ELSEIF ( MPI_COMPILER )
            # The mpi compiler should take care of everything
        ELSE()
            # Perform the default search for MPI
            INCLUDE ( FindMPI )
            IF ( NOT MPI_FOUND )
                MESSAGE( FATAL_ERROR "Did not find MPI" )
            ENDIF ()
            INCLUDE_DIRECTORIES ( ${MPI_INCLUDE_PATH} )
            SET ( MPI_INCLUDE ${MPI_INCLUDE_PATH} )
        ENDIF()
        # Check if we need to use MPI for serial tests
        CHECK_ENABLE_FLAG( USE_EXT_MPI_FOR_SERIAL_TESTS 0 )
        # Set defaults if they have not been set
        IF ( NOT MPIEXEC )
            SET( MPIEXEC mpirun )
        ENDIF()
        IF ( NOT MPIEXEC_NUMPROC_FLAG )
            SET( MPIEXEC_NUMPROC_FLAG "-np" )
        ENDIF()
        # Set the definitions
        ADD_DEFINITIONS( "-D USE_EXT_MPI" )  
        MESSAGE( "Using MPI" )
        MESSAGE( "  MPIEXEC = ${MPIEXEC}" )
        MESSAGE( "  MPIEXEC_NUMPROC_FLAG = ${MPIEXEC_NUMPROC_FLAG}" )
        MESSAGE( "  MPI_LINK_FLAGS = ${MPI_LINK_FLAGS}" )
        MESSAGE( "  MPI_LIBRARIES = ${MPI_LIBRARIES}" )
        MESSAGE( "  TEST_MAX_PROCS = ${TEST_MAX_PROCS}" )
    ELSE()
        SET( USE_EXT_MPI_FOR_SERIAL_TESTS 0 )
        SET( MPIEXEC "" )
        SET( MPIEXEC_NUMPROC_FLAG "" )
        SET( MPI_INCLUDE "" )
        SET( MPI_LINK_FLAGS "" )
        SET( MPI_LIBRARIES "" )
        MESSAGE( "Not using MPI, all parallel tests will be disabled" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure the libmesh libraries
MACRO ( CONFIGURE_LIBMESH )
    # Determine if we want to use libmesh
    CHECK_ENABLE_FLAG(USE_EXT_LIBMESH 1 )
    IF ( USE_EXT_LIBMESH )
        # Check if we specified the libmesh directory
        IF ( LIBMESH_DIRECTORY )
            LIBMESH_SET_INCLUDES( ${LIBMESH_DIRECTORY} )
            LIBMESH_SET_LIBRARIES( ${LIBMESH_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${LIBMESH_INCLUDE} )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for libmesh is not supported.  Use -D LIBMESH_DIRECTORY=" )
        ENDIF()
        MESSAGE( "Using libmesh" )
        MESSAGE( "   ${LIBMESH_LIBS}" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure NEK
MACRO ( CONFIGURE_NEK )
    # Determine if we want to use NEK
    CHECK_ENABLE_FLAG( USE_EXT_NEK 0 )
    IF ( USE_EXT_NEK )
        # Check if we specified the NEK directory
        IF ( NEK_DIRECTORY )
            VERIFY_PATH ( ${NEK_DIRECTORY} )
            # Include the NEK directories
            IF ( NOT NEK_INCLUDE )
                SET ( NEK_INCLUDE ${NEK_DIRECTORY} )
            ENDIF()
            # Find the NEK libaries
            IF ( NOT NEK_PATH_LIB )
                SET ( NEK_PATH_LIB ${NEK_DIRECTORY} )
            ENDIF()
            VERIFY_PATH ( ${NEK_PATH_LIB} )
            FIND_LIBRARY ( NEK_LIB     NAMES NEK5000      PATHS ${NEK_PATH_LIB}          NO_DEFAULT_PATH )
            IF ( NOT NEK_LIB )
                MESSAGE( FATAL_ERROR "Nek5000 library (NEK5000) not found in ${NEK_PATH_LIB}" )
            ENDIF ()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for NEK is not supported.  Use -D NEK_DIRECTORY=" )
        ENDIF()
        CHECK_ENABLE_FLAG( NOTIMER  0 )
        CHECK_ENABLE_FLAG( MPITIMER 0 )
        CHECK_ENABLE_FLAG( MPIIO    0 )
        CHECK_ENABLE_FLAG( BG       0 )
        CHECK_ENABLE_FLAG( K10_MXM  0 )
        CHECK_ENABLE_FLAG( CVODE    0 )
        CHECK_ENABLE_FLAG( NEKNEK   0 )
        CHECK_ENABLE_FLAG( MOAB     1 )
        IF ( NOT USE_EXT_MOAB ) 
            MESSAGE( FATAL_ERROR "Within AMP, MOAB is required to use Nek5000." )
        ENDIF()
        # Add the libraries in the appropriate order
        INCLUDE_DIRECTORIES ( ${NEK_INCLUDE} )
        SET ( NEK_LIBS
            ${NEK_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_NEK" )  
        MESSAGE( "Using NEK" )
        MESSAGE( "   ${NEK_LIBS}" )
        SET ( CURPACKAGE "nek" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure DENDRO
MACRO ( CONFIGURE_DENDRO )
    # Determine if we want to use 
    CHECK_ENABLE_FLAG( USE_EXT_DENDRO 0 )
    IF ( USE_EXT_DENDRO )
        IF ( DENDRO_DIRECTORY )
            VERIFY_PATH ( ${DENDRO_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${DENDRO_DIRECTORY}/include )
            SET ( DENDRO_INCLUDE ${DENDRO_DIRECTORY}/include )
            FIND_LIBRARY ( DENDRO_BIN_LIB   NAMES BinOps PATHS ${DENDRO_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( DENDRO_OCT_LIB   NAMES Oct    PATHS ${DENDRO_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( DENDRO_PAR_LIB   NAMES Par    PATHS ${DENDRO_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( DENDRO_POINT_LIB NAMES Point  PATHS ${DENDRO_DIRECTORY}/lib  NO_DEFAULT_PATH )
            FIND_LIBRARY ( DENDRO_TEST_LIB  NAMES Test   PATHS ${DENDRO_DIRECTORY}/lib  NO_DEFAULT_PATH )
            IF ( (NOT DENDRO_BIN_LIB) OR (NOT DENDRO_OCT_LIB) OR (NOT DENDRO_PAR_LIB) OR
                (NOT DENDRO_POINT_LIB) OR (NOT DENDRO_TEST_LIB) )
                MESSAGE( ${DENDRO_BIN_LIB} )
                MESSAGE( ${DENDRO_OCT_LIB} )
                MESSAGE( ${DENDRO_PAR_LIB} )
                MESSAGE( ${DENDRO_POINT_LIB} )
                MESSAGE( ${DENDRO_TEST_LIB} )
                MESSAGE( FATAL_ERROR "DENDRO libraries not found in ${DENDRO_DIRECTORY}/lib" )
            ENDIF ()
            # Add the libraries in the appropriate order
            SET ( DENDRO_LIBS
                ${DENDRO_OCT_LIB}
                ${DENDRO_PAR_LIB}
                ${DENDRO_POINT_LIB}
                ${DENDRO_TEST_LIB}
                ${DENDRO_BIN_LIB}
             )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for DENDRO is not supported.  Use -D DENDRO_DIRECTORY=" )
        ENDIF()
        ADD_DEFINITIONS ( "-D USE_EXT_DENDRO" )  
        MESSAGE( "Using DENDRO" )
        MESSAGE( "   ${DENDRO_LIBS}" )
    ENDIF()
ENDMACRO()


# Macro to find and configure MOAB
MACRO ( CONFIGURE_MOAB )
    # Determine if we want to use MOAB
    CHECK_ENABLE_FLAG( USE_EXT_MOAB 0 )
    IF ( USE_EXT_MOAB )
        # Check if we specified the MOAB directory
        IF ( MOAB_DIRECTORY )
            VERIFY_PATH ( ${MOAB_DIRECTORY} )
            # Include the MOAB directories
            SET ( MOAB_INCLUDE ${MOAB_DIRECTORY}/include )
            SET ( IMESH_INCLUDE ${MOAB_DIRECTORY}/lib )
            # Find the MOAB libaries
            SET ( MOAB_PATH_LIB ${MOAB_DIRECTORY}/lib )
            VERIFY_PATH ( ${MOAB_PATH_LIB} )
            FIND_LIBRARY ( MOAB_MESH_LIB     NAMES MOAB      PATHS ${MOAB_PATH_LIB}          NO_DEFAULT_PATH )
            FIND_LIBRARY ( MOAB_iMESH_LIB    NAMES iMesh     PATHS ${MOAB_PATH_LIB}          NO_DEFAULT_PATH )
            FIND_LIBRARY ( MOAB_COUPLER_LIB  NAMES mbcoupler PATHS ${MOAB_PATH_LIB}          NO_DEFAULT_PATH )
            IF ( NOT MOAB_MESH_LIB )
                MESSAGE( FATAL_ERROR "MOAB library (MOAB) not found in ${MOAB_PATH_LIB}" )
            ENDIF ()
            IF ( NOT MOAB_iMESH_LIB )
                MESSAGE( FATAL_ERROR "iMesh library ${MOAB_iMESH_LIB}  not found in ${MOAB_PATH_LIB}" )
            ENDIF ()
            IF ( NOT MOAB_COUPLER_LIB )
                MESSAGE( FATAL_ERROR "MBCoupler library ${MOAB_COUPLER_LIB}  not found in ${MOAB_PATH_LIB}" )
            ENDIF ()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for MOAB is not supported.  Use -D MOAB_DIRECTORY=" )
        ENDIF()
        # Check if we specified the cgm directory
        IF ( CGM_DIRECTORY )
            VERIFY_PATH ( ${CGM_DIRECTORY} )
            # Include the CGM directories
            SET ( MOAB_INCLUDE ${MOAB_INCLUDE} ${CGM_DIRECTORY}/include )
            # Find the CGM libaries
            SET ( CGM_PATH_LIB ${CGM_DIRECTORY}/lib )
            VERIFY_PATH ( ${CGM_PATH_LIB} )
            FIND_LIBRARY ( MOAB_CGM_LIB     NAMES cgm      PATHS ${CGM_PATH_LIB}        NO_DEFAULT_PATH )
            FIND_LIBRARY ( MOAB_iGEOM_LIB   NAMES iGeom    PATHS ${CGM_PATH_LIB}        NO_DEFAULT_PATH )
            IF ( NOT MOAB_CGM_LIB )
                MESSAGE( FATAL_ERROR "CGM library ${MOAB_CGM_LIB}  not found in ${CGM_PATH_LIB}" )
            ENDIF ()
            IF ( NOT MOAB_iGEOM_LIB )
                MESSAGE( FATAL_ERROR "iGEOM library ${MOAB_iGEOM_LIB}  not found in ${CGM_PATH_LIB}" )
            ENDIF ()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for cgm is not supported.  Use -D CGM_DIRECTORY=" )
        ENDIF()
        # Check if we specified the Cubit directory
        IF ( CUBIT_DIRECTORY )
            VERIFY_PATH ( ${CUBIT_DIRECTORY} )
            # Include the CUBIT directories
            # SET ( MOAB_INCLUDE ${MOAB_INCLUDE} ${CUBIT_DIRECTORY}/include )
            # Find the CGM libaries
            SET ( CUBIT_PATH_LIB ${CUBIT_DIRECTORY} )
            VERIFY_PATH ( ${CGM_PATH_LIB} )
            FIND_LIBRARY ( MOAB_CUBIT_LIB     NAMES cubiti19      PATHS ${CUBIT_PATH_LIB}        NO_DEFAULT_PATH )
            IF ( NOT MOAB_CUBIT_LIB )
                MESSAGE( FATAL_ERROR "CUBIT librarys not found in ${CUBIT_PATH_LIB}" )
            ENDIF ()
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for cubit is not supported.  Use -D CUBIT_DIRECTORY=" )
        ENDIF()
        # Add the libraries in the appropriate order
        INCLUDE_DIRECTORIES ( ${MOAB_INCLUDE} )
        SET ( MOAB_LIBS
            ${MOAB_COUPLER_LIB}
            ${MOAB_iMESH_LIB}
            ${MOAB_MESH_LIB}
            ${MOAB_CGM_LIB}
            ${MOAB_iGEOM_LIB}
            ${MOAB_CUBIT_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_MOAB" )  
        MESSAGE( "Using MOAB" )
        MESSAGE( "   ${MOAB_LIBS}" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure the sundials libraries
MACRO ( CONFIGURE_SUNDIALS_LIBRARIES )
    # Determine if we want to use sundials
    CHECK_ENABLE_FLAG(USE_EXT_SUNDIALS 1 )
    IF ( USE_EXT_SUNDIALS )
        # Check if we specified the libmesh directory
        IF ( SUNDIALS_DIRECTORY )
            SUNDIALS_SET_INCLUDES( ${SUNDIALS_DIRECTORY} )
            SUNDIALS_SET_LIBRARIES( ${SUNDIALS_DIRECTORY} )
            INCLUDE_DIRECTORIES ( ${SUNDIALS_INCLUDE} )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for sundials is not supported.  Use -D SUNDIALS_DIRECTORY=" )
        ENDIF()
        MESSAGE( "Using sundials" )
        MESSAGE( "   ${SUNDIALS_LIBS}" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure the hypre libraries
MACRO ( CONFIGURE_HYPRE_LIBRARIES )
    # Determine if we want to use silo
    CHECK_ENABLE_FLAG( USE_EXT_HYPRE 1 )
    IF ( USE_EXT_HYPRE )
        # Check if we specified the hypre directory
        IF ( HYPRE_DIRECTORY )
            VERIFY_PATH ( ${HYPRE_DIRECTORY} )
            SET ( HYPRE_LIB_DIRECTORY ${HYPRE_DIRECTORY}/lib )
            FIND_LIBRARY ( HYPRE_LIB         NAMES HYPRE                PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_DM_LIB      NAMES HYPRE_DistributedMatrix  PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_DMPS_LIB    NAMES HYPRE_DistributedMatrixPilutSolver  PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            # FIND_LIBRARY ( HYPRE_EUCLID_LIB  NAMES HYPRE_Euclid  PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_IJMV_LIB    NAMES HYPRE_IJ_mv          PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_KRYLOV_LIB  NAMES HYPRE_krylov         PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            # FIND_LIBRARY ( HYPRE_LSI_LIB     NAMES HYPRE_LSI  PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_MATMAT_LIB  NAMES HYPRE_MatrixMatrix   PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_MULTIV_LIB  NAMES HYPRE_multivector    PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_PARAS_LIB   NAMES HYPRE_ParaSails      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_PBMV_LIB    NAMES HYPRE_parcsr_block_mv  PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_PLS_LIB     NAMES HYPRE_parcsr_ls      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_PMV_LIB     NAMES HYPRE_parcsr_mv      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_SEQMV_LIB   NAMES HYPRE_seq_mv         PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_SSLS_LIB    NAMES HYPRE_sstruct_ls     PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_SSMV_LIB    NAMES HYPRE_sstruct_mv     PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_SLS_LIB     NAMES HYPRE_struct_ls      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_SMV_LIB     NAMES HYPRE_struct_mv      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
            FIND_LIBRARY ( HYPRE_UTIL_LIB    NAMES HYPRE_utilities      PATHS ${HYPRE_LIB_DIRECTORY}  NO_DEFAULT_PATH )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for hypre is not yet supported.  Use -D HYPRE_DIRECTORY=" )
        ENDIF()
        # Add the libraries in the appropriate order
        SET ( HYPRE_LIBS
            ${HYPRE_DM_LIB}
            ${HYPRE_DMPS_LIB}
            # ${HYPRE_EUCLID_LIB}
            ${HYPRE_IJMV_LIB}
            ${HYPRE_KRYLOV_LIB}
            # ${HYPRE_LSI_LIB}
            ${HYPRE_MATMAT_LIB}
            ${HYPRE_MULTIV_LIB}
            ${HYPRE_PARAS_LIB}
            ${HYPRE_PBMV_LIB}
            ${HYPRE_PLS_LIB}
            ${HYPRE_PMV_LIB}
            ${HYPRE_SEQMV_LIB}
            ${HYPRE_SSLS_LIB}
            ${HYPRE_SSMV_LIB}
            ${HYPRE_SLS_LIB}
            ${HYPRE_SMV_LIB}
            ${HYPRE_UTIL_LIB}
            ${HYPRE_LIB}
        )
        ADD_DEFINITIONS ( "-D USE_EXT_HYPRE" )  
        MESSAGE( "Using hypre" )
    ENDIF()
ENDMACRO ()


# Macro to find and configure the petsc libraries
MACRO ( CONFIGURE_PETSC_LIBRARIES )
    # Determine if we want to use petsc
    CHECK_ENABLE_FLAG(USE_EXT_PETSC 1 )
    IF ( USE_EXT_PETSC )
        # Check if we specified the petsc directory
        IF ( PETSC_DIRECTORY )
            VERIFY_PATH ( ${PETSC_DIRECTORY} )
            VERIFY_PATH ( ${PETSC_DIRECTORY}/include )
        ELSE()
            MESSAGE( FATAL_ERROR "Default search for petsc is not yet supported.  Use -D PETSC_DIRECTORY=" )
        ENDIF()
        # Get the petsc version
        PETSC_GET_VERSION( ${PETSC_DIRECTORY}/include )
        MESSAGE("Found PETSc version ${PETSC_VERSION}")
        # Add the petsc include folders and definitiosn
        SET ( PETSC_INCLUDE ${PETSC_INCLUDE} ${PETSC_DIRECTORY}/include )
        SET ( PETSC_INCLUDE ${PETSC_INCLUDE} ${PETSC_DIRECTORY}/${PETSC_ARCH}/include )
        SET ( PETSC_LIB_DIRECTORY ${PETSC_DIRECTORY}/${PETSC_ARCH}/lib )
        INCLUDE_DIRECTORIES ( ${PETSC_INCLUDE} )
        # Find the petsc libraries
        PETSC_SET_LIBRARIES( ${PETSC_LIB_DIRECTORY} )
        MESSAGE( "Using petsc" )
        MESSAGE( "   ${PETSC_LIBS}" )
    ENDIF()
ENDMACRO ()


# Macro to configure system-specific libraries and flags
MACRO ( CONFIGURE_SYSTEM )
    SET_COMPILER()
    # Remove extra library links
    CHECK_ENABLE_FLAG( USE_STATIC 0 )
    IF ( USE_STATIC )
        SET_STATIC_FLAGS()
    ENDIF()
    # Add system dependent flags
    MESSAGE("System is: ${CMAKE_SYSTEM_NAME}")
    IF ( ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
        # Windows specific system libraries
        SET( SYSTEM_PATHS "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.0A/Lib/x64" 
                          "C:/Program Files (x86)/Microsoft Visual Studio 8/VC/PlatformSDK/Lib/AMD64" )
        FIND_LIBRARY ( PSAPI_LIB    NAMES Psapi    PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        FIND_LIBRARY ( DBGHELP_LIB  NAMES DbgHelp  PATHS ${SYSTEM_PATHS}  NO_DEFAULT_PATH )
        SET( SYSTEM_LIBS ${PSAPI_LIB} ${DBGHELP_LIB} )
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux" )
        # Linux specific system libraries
        SET( SYSTEM_LIBS "-ldl" )
        CONFIGURE_ZLIB()
        IF ( NOT USE_STATIC )
            SET( SYSTEM_LIBS "${SYSTEM_LIBS} -rdynamic" )   # Needed for backtrace to print function names
        ENDIF()
        IF ( USING_GCC )
            SET( SYSTEM_LIBS "${SYSTEM_LIBS} -lgfortran" )   # Needed for backtrace to print function names
	ENDIF()
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin" )
        # Max specific system libraries
        SET( SYSTEM_LIBS "-ldl" )
        CONFIGURE_ZLIB()
        IF ( USING_GCC )
            SET( SYSTEM_LIBS "${SYSTEM_LIBS} -lgfortran" )
	ENDIF()
    ELSEIF( ${CMAKE_SYSTEM_NAME} STREQUAL "Generic" )
        # Generic system libraries
    ELSE()
        MESSAGE( FATAL_ERROR "OS not detected" )
    ENDIF()
        MESSAGE("System libs: ${SYSTEM_LIBS}")
ENDMACRO ()


# Macro to configure AMP-specific options
MACRO ( CONFIGURE_AMP )
    # Add the AMP install directory
    INCLUDE_DIRECTORIES ( ${AMP_INSTALL_DIR}/include )
    # Set the data directory for AMP (needed to find the meshes)
    IF ( AMP_DATA )
        VERIFY_PATH ( ${AMP_DATA} )
    ELSEIF ( NOT ONLY_BUILD_DOCS )
        MESSAGE( FATAL_ERROR "AMP_DATA must be set" )
    ENDIF()
    # Set the maximum number of processors for the tests
    IF ( NOT TEST_MAX_PROCS )
        SET( TEST_MAX_PROCS 32 )
    ENDIF()
    # Remove extra library links
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LINKER_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    # Check the user configure flags
    CHECK_ENABLE_FLAG( USE_AMP_UTILS 1 )
    CHECK_ENABLE_FLAG( USE_AMP_MESH 1 )
    CHECK_ENABLE_FLAG( USE_AMP_DISCRETIZATION 1 )
    CHECK_ENABLE_FLAG( USE_AMP_VECTORS 1 )
    CHECK_ENABLE_FLAG( USE_AMP_MATRICES 1 )
    CHECK_ENABLE_FLAG( USE_AMP_MATERIALS 1 )
    CHECK_ENABLE_FLAG( USE_AMP_OPERATORS 1 )
    CHECK_ENABLE_FLAG( USE_AMP_SOLVERS 1 )
    CHECK_ENABLE_FLAG( USE_AMP_TIME_INTEGRATORS 1 )
    # Check and disable packages based on dependencies
    IF ( NOT BUILD_ONLY_DOCS )
        # Check if we are using utils
        IF ( NOT USE_AMP_UTILS )
            MESSAGE( FATAL_ERROR "AMP Utils must be used" )
        ENDIF()
        # Check if we are using ampmesh
        IF ( NOT USE_AMP_MESH )
            MESSAGE( "Disabling AMP Mesh" )
        ENDIF()
        # Check if we are using discretization
        IF ( NOT USE_AMP_MESH )
            SET ( USE_AMP_DISCRETIZATION 0 )
        ENDIF()
        IF ( NOT USE_AMP_DISCRETIZATION )
            MESSAGE( "Disabling AMP Descritization" )
        ENDIF()
        # Check if we are using vectors
        IF ( NOT USE_AMP_VECTORS )
            MESSAGE( "Disabling AMP Vectors" )
        ENDIF()
        # Check if we are using matrices
        IF ( NOT USE_AMP_MATRICES )
            MESSAGE( "Disabling AMP Matrices" )
        ENDIF()
        # Check if we are using materials
        IF ( NOT USE_AMP_VECTORS )
            SET ( USE_AMP_MATERIALS 0 )
        ENDIF()
        IF ( NOT USE_AMP_MATRICES )
            MESSAGE( "Disabling AMP Materials" )
        ENDIF()
        # Check if we are using operators
        IF ( (NOT USE_AMP_MESH) OR (NOT USE_AMP_VECTORS) OR (NOT USE_AMP_MATRICES) )
            SET ( USE_AMP_OPERATORS 0 )
        ENDIF()
        IF ( NOT USE_AMP_OPERATORS )
            MESSAGE( "Disabling AMP Operators" )
        ENDIF()
        # Check if we are using solvers
        IF ( (NOT USE_AMP_OPERATORS) OR (NOT USE_EXT_LIBMESH) )
            SET ( USE_AMP_SOLVERS 0 )
        ENDIF()
        IF ( NOT USE_AMP_SOLVERS )
            MESSAGE( "Disabling AMP Solvers" )
        ENDIF()
        # Check if we are using time_integrators
        IF ( (NOT USE_AMP_SOLVERS) )
            SET ( USE_AMP_TIME_INTEGRATORS 0 )
        ENDIF()
        IF ( NOT USE_AMP_TIME_INTEGRATORS )
            MESSAGE( "Disabling AMP Time Integrators" )
        ENDIF()
    ENDIF()
    # Set which packages we enabled
    SET( ${PROJECT_NAME}_ENABLE_AMP_UTILITIES       OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_MESH            OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_DISCRETIZATION  OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_VECTORS         OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_MATRICES        OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_MATERIALS       OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_OPERATORS       OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_TIME_INTEGRATORS OFF )
    SET( ${PROJECT_NAME}_ENABLE_AMP_SOLVERS         OFF )
    IF ( USE_AMP_UTILS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_UTILITIES ON )
        SET( AMP_UTILITIES_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_MESH )
        SET( ${PROJECT_NAME}_ENABLE_AMP_MESH ON )
        SET( AMP_MESH_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_DISCRETIZATION )
        SET( ${PROJECT_NAME}_ENABLE_AMP_DISCRETIZATION ON )
        SET( AMP_DISCRETIZATION_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_VECTORS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_VECTORS ON )
        SET( AMP_VECTORS_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_MATRICES )
        SET( ${PROJECT_NAME}_ENABLE_AMP_MATRICES ON )
        SET( AMP_MATRICES_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_MATERIALS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_MATERIALS ON )
        SET( AMP_MATERIALS_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_OPERATORS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_OPERATORS ON )
        SET( AMP_OPERATORS_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_TIME_INTEGRATORS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_TIME_INTEGRATORS ON )
        SET( AMP_TIME_INTEGRATORS_ENABLE_TESTS ON )
    ENDIF()
    IF ( USE_AMP_SOLVERS )
        SET( ${PROJECT_NAME}_ENABLE_AMP_SOLVERS ON )
        SET( AMP_SOLVERS_ENABLE_TESTS ON )
    ENDIF()
    GLOBAL_SET( USE_AMP_UTILS "" )
    GLOBAL_SET( USE_AMP_MESH "" )
    GLOBAL_SET( USE_AMP_DISCRETIZATION "" )
    GLOBAL_SET( USE_AMP_VECTORS "" )
    GLOBAL_SET( USE_AMP_MATRICES "" )
    GLOBAL_SET( USE_AMP_MATERIALS "" )
    GLOBAL_SET( USE_AMP_OPERATORS "" )
    GLOBAL_SET( USE_AMP_TIME_INTEGRATORS "" )
    GLOBAL_SET( USE_AMP_SOLVERS "" )
    GLOBAL_SET ( AMP_DOC_DIRS " ")
    IF ( USE_EXT_NEK )
        SET ( AMP_LIBS ${AMP_LIBS} "nek" )
        ADD_DEFINITIONS ( -D USE_EXT_NEK )  
    ENDIF()
    # Add documentation folders
    IF ( USE_AMP_TIME_INTEGRATORS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/time_integrators\"" )
    ENDIF()
    IF ( USE_AMP_SOLVERS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/solvers\"" )
    ENDIF()
    IF ( USE_AMP_OPERATORS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/operators\"" )
    ENDIF()
    IF ( USE_AMP_MATERIALS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/materials\"" )
    ENDIF()
    IF ( USE_AMP_MATRICES )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/matrices\"" )
    ENDIF()
    IF ( USE_AMP_VECTORS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/vectors\"" )
    ENDIF()
    IF ( USE_AMP_DISCRETIZATION )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/discretization\"" )
    ENDIF()
    IF ( USE_AMP_MESH )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/ampmesh\"" )
    ENDIF()
    IF ( USE_AMP_UTILS )
        SET ( AMP_DOC_DIRS "${AMP_DOC_DIRS}  \"${AMP_SOURCE_DIR}/src/utils\"" )
    ENDIF()
    # Set the flags
    SET_AMP_PACKAGE_FLAGS()
ENDMACRO ()



