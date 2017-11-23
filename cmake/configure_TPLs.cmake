# Macro to configure the TPL builder
FUNCTION ( CONFIGURE_TPL_BUILDER )
    # Set the build type
    SET( TPL_CMAKE "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}" )
    # Set the compiler / compile flags
    IF ( C_COMPILER )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DC_COMPILER=${C_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCFLAGS=${CFLAGS}" )
    ELSE()
        SET( TPL_CMAKE ${TPL_CMAKE} "-DC_COMPILER=${CMAKE_C_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCFLAGS=${CMAKE_C_FLAGS}" )
    ENDIF()
    IF ( CXX_COMPILER )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCXX_COMPILER=${CXX_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCXXFLAGS=${CXXFLAGS}" )
    ELSE()
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DCXXFLAGS=${CMAKE_CXX_FLAGS}" )
    ENDIF()
    IF ( Fortran_COMPILER )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DFortran_COMPILER=${Fortran_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DFFLAGS=${FFLAGS}" )
    ELSE()
        SET( TPL_CMAKE ${TPL_CMAKE} "-DFortran_COMPILER=${CMAKE_Fortran_COMPILER}" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DFFLAGS=${CMAKE_Fortran_FLAGS}" )
    ENDIF()
    IF ( ENABLE_STATIC )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_STATIC:BOOL=ON" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_SHARED:BOOL=OFF" )
    ELSEIF()
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_STATIC:BOOL=OFF" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_SHARED:BOOL=ON" )
    ELSE()
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_STATIC:BOOL=ON" )
        SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_SHARED:BOOL=OFF" )
    ENDIF()
    SET( TPL_CMAKE ${TPL_CMAKE} "-DCXX_STD=${CXX_STD}" )
    SET( TPL_CMAKE ${TPL_CMAKE} "-DLDFLAGS=${LDFLAGS}" )
    # Set the install path
    SET( TPL_CMAKE ${TPL_CMAKE} "-DINSTALL_DIR:PATH=${${PROJ}_INSTALL_DIR}/TPLs" )
    SET( TPL_CMAKE ${TPL_CMAKE} "-DENABLE_TESTS:BOOL=OFF" )
    # Set all flags that start with TPL_
    GET_CMAKE_PROPERTY( variableNames VARIABLES )
    FOREACH ( var ${variableNames} )
        STRING( REGEX REPLACE "^TPL_" "" var2 "${var}" )
        IF ( ${var} STREQUAL TPL_URL )
            # Special variable
        ELSEIF ( ${var} STREQUAL TPL_LIST )
            STRING( REPLACE ";" "," TPL_LIST2 "${TPL_LIST}" )
            SET( TPL_CMAKE ${TPL_CMAKE} "-DTPL_LIST:STRING=${TPL_LIST2}" )
        ELSEIF ( ${var} STREQUAL TPL_${var2} )
            SET( TPL_CMAKE ${TPL_CMAKE} "-D${var2}=${${var}}" )
        ENDIF()
    ENDFOREACH()
    # Download the TPL builder
    MESSAGE( STATUS "Downloading TPL builder" )
    IF ( NOT EXISTS tpl-builder )
        EXECUTE_PROCESS(
            COMMAND hg clone "${TPL_URL}"
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
            OUTPUT_QUIET
        )
    ENDIF()
    EXECUTE_PROCESS(
        COMMAND hg pull -u
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tpl_directory"
        OUTPUT_QUIET
    )
    # Configure the TPL builder
    MESSAGE( STATUS "Configuring TPL builder" )
    FILE( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tpl-build" )
    EXECUTE_PROCESS(
        COMMAND ${CMAKE_COMMAND} ${TPL_CMAKE} "${CMAKE_BINARY_DIR}/tpl-builder"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tpl-build"
    )
    # Build the TPLs
    MESSAGE( STATUS "Building TPLs" )
    EXECUTE_PROCESS(
        COMMAND ${CMAKE_MAKE_PROGRAM}
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/tpl-build"
    )
    SET( TPL_DIRECTORY "${${PROJ}_INSTALL_DIR}/TPLs" PARENT_SCOPE )
ENDFUNCTION()
