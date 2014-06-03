INCLUDE(CheckCSourceCompiles)


MACRO(GLOBAL_SET VARNAME)
  SET(${VARNAME} ${ARGN} CACHE INTERNAL "")
ENDMACRO()


# CMake assert
MACRO(ASSERT test comment)
    IF (NOT ${test})
        MESSSAGE(FATAL_ERROR "Assertion failed: ${comment}")
    ENDIF(NOT ${test})
ENDMACRO(ASSERT)


# Macro to convert a m4 file
# This command converts a file of the format "global_path/file.fm4"
# and convertes it to file.f90.  It also requires the path.  
MACRO (CONVERT_M4_FORTRAN IN LOCAL_PATH OUT_PATH )
    STRING(REGEX REPLACE ${LOCAL_PATH} "" OUT ${IN} )
    STRING(REGEX REPLACE "/" "" OUT ${OUT} )
    STRING(REGEX REPLACE ".fm4" ".f90" OUT ${CMAKE_CURRENT_BINARY_DIR}/${OUT_PATH}/${OUT} )
    CONFIGURE_FILE ( ${IN} ${IN} COPYONLY )
    add_custom_command(
        OUTPUT ${OUT}
        COMMAND m4 --prefix-builtins -I${LOCAL_PATH} ${M4_OPTIONS} ${IN} > ${OUT}
        DEPENDS ${IN}
    )
    set_source_files_properties(${OUT} PROPERTIES GENERATED true)
    SET( SOURCES ${SOURCES} "${OUT}" )
ENDMACRO ()


# Add a package to the test dependency list
MACRO (ADD_PACKAGE_TO_TEST_DEP_LIST PACKAGE)
    IF ( TEST_DEP_LIST )
        SET( TEST_DEP_LIST ${PACKAGE} ${TEST_DEP_LIST} )
    ELSE()
        SET( TEST_DEP_LIST ${PACKAGE} )
    ENDIF()
ENDMACRO ()


# Add a package to the AMP library
MACRO( ADD_AMP_LIBRARY PACKAGE )
    ADD_PACKAGE_TO_TEST_DEP_LIST( ${PACKAGE} )
    #INCLUDE_DIRECTORIES( ${AMP_INSTALL_DIR}/include/${PACKAGE} )
    ADD_SUBDIRECTORY( ${PACKAGE} )
ENDMACRO()


# Add an executable
MACRO( ADD_AMP_EXECUTABLE PACKAGE )
    ADD_SUBDIRECTORY( ${PACKAGE} )
ENDMACRO()


# Initialize a package
MACRO (BEGIN_PACKAGE_CONFIG PACKAGE)
    SET( HEADERS "" )
    SET( CXXSOURCES "" )
    SET( CSOURCES "" )
    SET( FSOURCES "" )
    SET( M4FSOURCES "" )
    SET( SOURCES "" )
    SET( CURPACKAGE ${PACKAGE} )
ENDMACRO ()


# Find the source files
MACRO (FIND_FILES)
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "*.h" "*.hh" "*.hpp" "*.I" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "*.cc" "*.cpp" "*.cxx" "*.C" )
    # Find the Fortran sources
    SET( T_FSOURCES "" )
    FILE( GLOB T_FSOURCES "*.f" "*.f90" )
    # Find the m4 fortran source (and convert)
    SET( T_M4FSOURCES "" )
    FILE( GLOB T_M4FSOURCES "*.fm4" )
    FOREACH( m4file ${T_M4FSOURCES} )
        CONVERT_M4_FORTRAN( ${m4file} ${CMAKE_CURRENT_SOURCE_DIR} "" )
    ENDFOREACH ()
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( FSOURCES ${FSOURCES} ${T_FSOURCES} )
    SET( M4FSOURCES ${M4FSOURCES} ${T_M4FSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_FSOURCES} ${T_M4FSOURCES} )
ENDMACRO()


# Find the source files
MACRO (FIND_FILES_PATH IN_PATH)
    # Find the C/C++ headers
    SET( T_HEADERS "" )
    FILE( GLOB T_HEADERS "${IN_PATH}/*.h" "${IN_PATH}/*.hh" "${IN_PATH}/*.hpp" "${IN_PATH}/*.I" )
    # Find the C sources
    SET( T_CSOURCES "" )
    FILE( GLOB T_CSOURCES "${IN_PATH}/*.c" )
    # Find the C++ sources
    SET( T_CXXSOURCES "" )
    FILE( GLOB T_CXXSOURCES "${IN_PATH}/*.cc" "${IN_PATH}/*.cpp" "${IN_PATH}/*.cxx" )
    # Find the Fortran sources
    SET( T_FSOURCES "" )
    FILE( GLOB T_FSOURCES "${IN_PATH}/*.f" "${IN_PATH}/*.f90" )
    # Find the m4 fortran source (and convert)
    SET( T_M4FSOURCES "" )
    FILE( GLOB T_M4FSOURCES "${IN_PATH}/*.fm4" )
    FOREACH (m4file ${T_M4FSOURCES})
        CONVERT_M4_FORTRAN( ${m4file} ${CMAKE_CURRENT_SOURCE_DIR}/${IN_PATH} "" )
    ENDFOREACH ()
    # Add all found files to the current lists
    SET( HEADERS ${HEADERS} ${T_HEADERS} )
    SET( CXXSOURCES ${CXXSOURCES} ${T_CXXSOURCES} )
    SET( CSOURCES ${CSOURCES} ${T_CSOURCES} )
    SET( FSOURCES ${FSOURCES} ${T_FSOURCES} )
    SET( SOURCES ${SOURCES} ${T_CXXSOURCES} ${T_CSOURCES} ${T_FSOURCES} )
ENDMACRO()


# Add a subdirectory
MACRO( ADD_PACKAGE_SUBDIRECTORY SUBDIR )
    FIND_FILES_PATH( ${SUBDIR} )
ENDMACRO()


# Install a package
MACRO( INSTALL_AMP_TARGET LIBNAME )
    # Find all files in the current directory
    FIND_FILES()
    # Copy the header files to the include path
    FILE( GLOB HFILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${HEADERS} )
    STRING(REGEX REPLACE "${AMP_SOURCE_DIR}/" "" COPY_TARGET "copy-AMP-${CMAKE_CURRENT_SOURCE_DIR}-include" )
    STRING(REGEX REPLACE "/" "-" COPY_TARGET ${COPY_TARGET} )
    IF(NOT TARGET ${COPY_TARGET})
        ADD_CUSTOM_TARGET ( ${COPY_TARGET} ALL )
        ADD_DEPENDENCIES ( copy-AMP-Data ${COPY_TARGET} )
    ENDIF()
    FOREACH( HFILE ${HFILES} )
        SET( SRC_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${HFILE}" )
        SET( DST_FILE "${AMP_INSTALL_DIR}/include/${CURPACKAGE}/${HFILE}" )
        ADD_CUSTOM_COMMAND(TARGET ${COPY_TARGET} 
            PRE_BUILD 
            COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_FILE}" "${DST_FILE}"
        )
        ADD_DEPENDENCIES( ${COPY_TARGET}  "${CMAKE_CURRENT_SOURCE_DIR}/${HFILE}" )
    ENDFOREACH()
    ADD_DEPENDENCIES( copy-AMP-include ${COPY_TARGET} )
    # Add the library and install the package
    IF ( NOT ONLY_BUILD_DOCS )
        ADD_LIBRARY( ${LIBNAME} ${SOURCES} )
        ADD_DEPENDENCIES ( ${LIBNAME} copy-AMP-include )
        INSTALL( TARGETS ${LIBNAME} DESTINATION ${AMP_INSTALL_DIR}/lib )
    ENDIF()
    INSTALL( FILES ${HEADERS} DESTINATION ${AMP_INSTALL_DIR}/include/${LIBNAME} )
ENDMACRO()


# Macro to verify that a variable has been set
MACRO( VERIFY_VARIABLE VARIABLE_NAME )
    IF ( NOT ${VARIABLE_NAME} )
        MESSAGE( FATAL_ERROR "PLease set: " ${VARIABLE_NAME} )
    ENDIF()
ENDMACRO()


# Macro to verify that a path has been set
MACRO( VERIFY_PATH PATH_NAME )
    IF ("${PATH_NAME}" STREQUAL "")
        MESSAGE ( FATAL_ERROR "Path is not set: ${PATH_NAME}" )
    ENDIF()
    IF ( NOT EXISTS ${PATH_NAME} )
        MESSAGE ( FATAL_ERROR "Path does not exist: ${PATH_NAME}" )
    ENDIF()
ENDMACRO()


# Macro to tell cmake to use static libraries
MACRO( SET_STATIC_FLAGS )
    # Remove extra library links
    set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
    set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
    set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
    set(CMAKE_SHARED_LINKER_FLAGS)
    set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
    set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
    # Add the static flag if necessary
    SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")    # Add static flag
    SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "-static")  # Add static flag
ENDMACRO()

# Macro to identify the compiler
MACRO( SET_COMPILER )
    # SET the C/C++ compiler
    IF( CMAKE_COMPILE_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX )
        SET( USING_GCC TRUE )
        MESSAGE("Using gcc")
    ELSEIF( MSVC OR MSVC_IDE OR MSVC60 OR MSVC70 OR MSVC71 OR MSVC80 OR CMAKE_COMPILER_2005 OR MSVC90 OR MSVC10 )
        IF( NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows" )
            MESSAGE( FATAL_ERROR "Using microsoft compilers on non-windows system?" )
        ENDIF()
        SET( USING_MICROSOFT TRUE )
        MESSAGE("Using Microsoft")
    ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "Intel") OR (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel") ) 
        SET(USING_ICC TRUE)
        MESSAGE("Using icc")
    ELSEIF( ${CMAKE_C_COMPILER_ID} MATCHES "PGI")
        SET(USING_PGCC TRUE)
        MESSAGE("Using pgCC")
    ELSEIF( (${CMAKE_C_COMPILER_ID} MATCHES "CRAY") OR (${CMAKE_C_COMPILER_ID} MATCHES "Cray") )
        SET(USING_CRAY TRUE)
        MESSAGE("Using Cray")
    ELSE()
        SET(USING_DEFAULT TRUE)
        MESSAGE("${CMAKE_C_COMPILER_ID}")
        MESSAGE("Unknown C/C++ compiler, default flags will be used")
    ENDIF()
    # SET the Fortran++ compiler
    IF( USE_FORTRAN )
        IF( CMAKE_COMPILE_IS_GFORTRAN OR (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU") )
            SET( USING_GFORTRAN TRUE )
            MESSAGE("Using gfortran")
        ELSEIF( (${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel") ) 
            SET(USING_IFORT TRUE)
            MESSAGE("Using ifort")
        ELSEIF( ${CMAKE_Fortran_COMPILER_ID} MATCHES "PGI")
            SET(USING_PGF90 TRUE)
            MESSAGE("Using pgf90")
        ELSE()
            SET(USING_DEFAULT TRUE)
            MESSAGE("${CMAKE_Fortran_COMPILER_ID}")
            MESSAGE("Unknown Fortran compiler, default flags will be used")
        ENDIF()
    ENDIF()
ENDMACRO()


# Macro to set the proper warnings
MACRO ( SET_WARNINGS )
  IF ( USING_GCC )
    # Add gcc specific compiler options
    #    -Wno-reorder:  warning: "" will be initialized after "" when initialized here
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} -Wall ") 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall ")
    # Disable warnings that I think are irrelavent (may need to be revisited)
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wno-reorder " )
    # Disable warnings that occur frequently, but should be fixed eventually
    SET(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS} -Wno-unused-variable" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wno-unused-variable" )
    # Add gcc specific flags
    SET(CMAKE_C_FLAGS " ${CMAKE_C_FLAGS}" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}" )
  ELSEIF ( USING_MICROSOFT )
    # Add Microsoft specifc compiler options
    SET(CMAKE_C_FLAGS     "${CMAKE_C_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D _SCL_SECURE_NO_WARNINGS /D _CRT_SECURE_NO_WARNINGS /D _ITERATOR_DEBUG_LEVEL=0" )
  ELSEIF ( USING_ICC )
    # Add Intel specifc compiler options
    #    111: statement is unreachable
    #         This occurs in LibMesh
    #    177: variable "" was declared but never referenced
    #    181: argument is incompatible with corresponding format string conversion
    #    304: access control not specified ("public" by default)
    #         This occurs in LibMesh
    #    383: value copied to temporary, reference to temporary used
    #         This is an irrelavent error
    #    444: destructor for base class "" is not virtual
    #         This can create memory leaks (and should be fixed)
    #         Unfortunatelly many of these come from LibMesh
    #    522: function "xxx" redeclared "inline" after being called
    #         We should fix this, but there are a lot of these
    #    593: variable "xxx" was set but never used
    #         This occurs commonly for error flags and internal variables that are helpful for debugging
    #    654: overloaded virtual function "" is only partially overridden in class " "
    #    869: parameter "xxx" was never referenced
    #         I believe this is bad practice and should be fixed, but it may require a broader discussion (it is built into the design of Operator)
    #    981: operands are evaluated in unspecified order
    #         This can occur when an implicit conversion take place in a function call 
    #   1011: missing return statement at end of non-void function
    #         This is bad practice
    #   1418: external function definition with no prior declaration
    #         This can happen if we don't include a header file (and maybe if there is an internal function?)
    #         Unfortunatelly many of these come from trilinos
    #   1419: external declaration in primary source file
    #         This occurs in a lot of the input processor, and needs to be revisited
    #   1572: floating-point equality and inequality comparisons are unreliable
    #         LibMesh warnings
    #   1599: declaration hides parameter 
    #         LibMesh warnings
    #   2259: non-pointer conversion from "int" to "unsigned char" may lose significant bits
    #         This is bad practice, use an explict coversion instead
    #         Unfortunatelly many of these come from LibMesh
    #   6843: A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall" )
    # Disable warnings that I think are irrelavent (may need to be revisited)
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -wd383 -wd593 -wd981" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -wd383 -wd593 -wd981" )
    # Disable warnings that occur due to other packages (it would be nice to disable them for certain header files only)
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -wd111 -wd304 -wd304 -wd444 -wd1418 -wd1572 -wd1599 -wd2259" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -wd111 -wd304 -wd304 -wd444 -wd1418 -wd1572 -wd1599 -wd2259" )
    SET(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -diag-disable 6843" )
    # Disable warnings that occur frequently, but should be fixed eventually
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -wd522 -wd869 -wd1419" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -wd522 -wd869 -wd1419" )
  ELSEIF ( USING_CRAY )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}")
  ELSEIF ( USING_PGCC )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS}")
  ELSEIF ( USING_DEFAULT )
    # Add default compiler options
    SET(CMAKE_C_FLAGS     " ${CMAKE_C_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -Wall")
  ENDIF ()
ENDMACRO ()


# Macro to add user compile flags
MACRO( ADD_USER_FLAGS )
    SET(CMAKE_C_FLAGS   " ${CMAKE_C_FLAGS} ${CFLAGS}" )
    SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${CXXFLAGS}" )
    SET(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} ${FFLAGS}" )
ENDMACRO()


# Macro to set the flags for debug mode
MACRO( SET_COMPILER_FLAGS )
    # Initilaize the compiler
    SET_COMPILER()
    # Set the default flags for each build type
    IF ( USING_MICROSOFT )
        SET(CMAKE_C_FLAGS_DEBUG       "-D_DEBUG /DEBUG /Od /EHsc /MDd /Zi" )
        SET(CMAKE_C_FLAGS_RELEASE     "/O2 /EHsc /MD"                      )
        SET(CMAKE_CXX_FLAGS_DEBUG     "-D_DEBUG /DEBUG /Od /EHsc /MDd /Zi" )
        SET(CMAKE_CXX_FLAGS_RELEASE   "/O2 /EHsc /MD"                      )
        SET(CMAKE_Fortran_FLAGS_DEBUG ""                                   )
        SET(CMAKE_Fortran_FLAGS_RELEASE ""                                 )
    ELSE()
        SET(CMAKE_C_FLAGS_DEBUG       "-g -D_DEBUG -O0" )
        SET(CMAKE_C_FLAGS_RELEASE     "-O2"             )
        SET(CMAKE_CXX_FLAGS_DEBUG     "-g -D_DEBUG -O0" )
        SET(CMAKE_CXX_FLAGS_RELEASE   "-O2"             )
        SET(CMAKE_Fortran_FLAGS_DEBUG "-g -O0"          )
        SET(CMAKE_Fortran_FLAGS_RELEASE "-O2"           )
    ENDIF()
    IF ( NOT DISABLE_GXX_DEBUG )
        SET(CMAKE_C_FLAGS_DEBUG   " ${CMAKE_C_FLAGS_DEBUG}   -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
        SET(CMAKE_CXX_FLAGS_DEBUG " ${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC" )
    ENDIF ()
    # Set the compiler flags to use
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR ${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
        SET(CMAKE_C_FLAGS       ${CMAKE_C_FLAGS_DEBUG}       )
        SET(CMAKE_CXX_FLAGS     ${CMAKE_CXX_FLAGS_DEBUG}     )
        SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG} )
    ELSEIF ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" OR ${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
        SET(CMAKE_C_FLAGS       ${CMAKE_C_FLAGS_RELEASE}       )
        SET(CMAKE_CXX_FLAGS     ${CMAKE_CXX_FLAGS_RELEASE}     )
        SET(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE} )
    ELSE()
        MESSAGE(FATAL_ERROR "Unknown value for CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
    ENDIF()
    # Add the user flags
    ADD_USER_FLAGS()
    # Set the warnings to use
    SET_WARNINGS()
ENDMACRO()


# Macro to copy data file at build time
MACRO(COPY_DATA_FILE SRC_FILE DST_FILE)
    STRING(REGEX REPLACE "${AMP_SOURCE_DIR}/" "" COPY_TARGET "copy-AMP-${CMAKE_CURRENT_SOURCE_DIR}" )
    STRING(REGEX REPLACE "/" "-" COPY_TARGET ${COPY_TARGET} )
    if(NOT TARGET ${COPY_TARGET})
        ADD_CUSTOM_TARGET ( ${COPY_TARGET} ALL )
        ADD_DEPENDENCIES ( copy-AMP-Data ${COPY_TARGET} )
    ENDIF()
    ADD_CUSTOM_COMMAND(TARGET ${COPY_TARGET} 
        PRE_BUILD 
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SRC_FILE}" "${DST_FILE}"
    )
ENDMACRO()


# Macro to copy a data file
MACRO ( COPY_TEST_DATA_FILE FILENAME ${ARGN} )
    SET( FILE_TO_COPY  ${CMAKE_CURRENT_SOURCE_DIR}/data/${FILENAME} )
    FOREACH( tmp ${ARGN})
        SET( FILE_TO_COPY  ${CMAKE_CURRENT_SOURCE_DIR}/${tmp}/${FILENAME} )
    ENDFOREACH()
    SET( DESTINATION_NAME ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
    IF ( EXISTS ${FILE_TO_COPY} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION_NAME} )
    ELSE()
        MESSAGE ( WARNING "Cannot find file: " ${FILE_TO_COPY} )
    ENDIF()
ENDMACRO ()


# Macro to copy a data file
MACRO ( COPY_EXAMPLE_DATA_FILE FILENAME )
    SET( FILE_TO_COPY  ${CMAKE_CURRENT_SOURCE_DIR}/data/${FILENAME} )
    SET( DESTINATION1 ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
    SET( DESTINATION2 ${EXAMPLE_INSTALL_DIR}/${FILENAME} )
    IF ( EXISTS ${FILE_TO_COPY} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION1} )
        COPY_DATA_FILE( ${FILE_TO_COPY} ${DESTINATION2} )
    ELSE()
        MESSAGE( WARNING "Cannot find file: " ${FILE_TO_COPY} )
    ENDIF()
ENDMACRO()


# Macro to copy a mesh file
MACRO (COPY_MESH_FILE MESHNAME)
    # Check the local data directory
    FILE ( GLOB MESHPATH ${CMAKE_CURRENT_SOURCE_DIR}/data/${MESHNAME} )
    # Check the AMP_DATA/meshes directory
    IF ( NOT MESHPATH )
        FILE ( GLOB MESHPATH ${AMP_DATA}/meshes/TestMeshes/${MESHNAME} )
    ENDIF ()
    IF ( NOT MESHPATH )
        FILE ( GLOB_RECURSE MESHPATH ${AMP_DATA}/meshes/*/${MESHNAME} )
    ENDIF ()
    # Check the AMP_DATA/vvu directory
    IF ( NOT MESHPATH )
        FILE ( GLOB MESHPATH ${AMP_DATA}/vvu/meshes/${MESHNAME} )
    ENDIF ()
    # We have either found the mesh or failed
    IF ( NOT MESHPATH )
        MESSAGE ( WARNING "Cannot find mesh: " ${MESHNAME} )
    ELSE ()
        STRING(REGEX REPLACE ${MESHNAME} "" MESHPATH ${MESHPATH} )
        COPY_DATA_FILE( ${MESHPATH}/${MESHNAME} ${CMAKE_CURRENT_BINARY_DIR}/${MESHNAME} )
    ENDIF ()
ENDMACRO()


# Macro to add the dependencies and libraries to an executable
MACRO ( ADD_AMP_EXE_DEP EXE )
    # Add the package dependencies
    IF ( AMP_TEST_LIB_EXISTS )
        ADD_DEPENDENCIES ( ${EXE} ${PACKAGE_TEST_LIB} )
        TARGET_LINK_LIBRARIES ( ${EXE} ${PACKAGE_TEST_LIB} )
    ENDIF()
    # Add the executable to the dependencies of check and build-test
    ADD_DEPENDENCIES( check ${EXE} )
    ADD_DEPENDENCIES( build-test ${EXE} )
    # Add the file copy targets to the dependency list
    ADD_DEPENDENCIES ( ${EXE} copy-AMP-Data )
    # Add the libraries
    TARGET_LINK_LIBRARIES ( ${EXE} ${AMP_LIBS} ${AMP_LIBS} )
    TARGET_LINK_LIBRARIES ( ${EXE} ${${PROJECT_NAME}_LIBRARIES} )
    # Add external libraries
    SET_TARGET_PROPERTIES(${EXE} PROPERTIES LINK_FLAGS "${LDFLAGS}")
    TARGET_LINK_LIBRARIES ( ${EXE} ${LIBMESH_LIBS} ${NEK_LIBS} ${MOAB_LIBS} ${DENDRO_LIBS} ${TRILINOS_LIBS} ${NETCDF_LIBS} ${PETSC_LIBS} ${X11_LIBS} ${SILO_LIBS} ${HDF5_LIBS} ${HYPRE_LIBS} )
    IF ( ${USE_EXT_SUNDIALS} )
        TARGET_LINK_LIBRARIES ( ${EXE} ${SUNDIALS_LIBS} )
    ENDIF  ()
    TARGET_LINK_LIBRARIES ( ${EXE} ${MPI_LINK_FLAGS} ${MPI_LIBRARIES} )
    TARGET_LINK_LIBRARIES ( ${EXE} ${LAPACK_LIBS} ${BLAS_LIBS} ${BLAS_LAPACK_LIBS} )
    TARGET_LINK_LIBRARIES ( ${EXE} ${COVERAGE_LIBS} ${LDLIBS} )
    TARGET_LINK_LIBRARIES ( ${EXE} ${SYSTEM_LIBS} )
ENDMACRO()


# Add a executable
MACRO( INSTALL_AMP_EXE EXE )
    SET( SOURCES ${EXE}.cc )
    ADD_EXECUTABLE( ${EXE} ${SOURCES} )
    ADD_AMP_EXE_DEP( ${EXE} )
ENDMACRO()


MACRO( ADD_FILES_TO_TEST_LIB FILENAMES )
    ADD_LIBRARY( ${PACKAGE_TEST_LIB} ${FILENAMES} )
    IF ( TEST_DEP_LIST )
        TARGET_LINK_LIBRARIES ( ${PACKAGE_TEST_LIB} ${TEST_DEP_LIST} )
    ENDIF()
    SET( AMP_TEST_LIB_EXISTS ${PACKAGE_TEST_LIB} )
ENDMACRO()


# Check if we want to keep the test
FUNCTION( KEEP_TEST RESULT )
    SET( ${RESULT} 1 PARENT_SCOPE )
    IF ( NOT ${PACKAGE_NAME}_ENABLE_TESTS )
        SET( ${RESULT} 0 PARENT_SCOPE )
    ENDIF()
    IF ( ONLY_BUILD_DOCS )
        SET( ${RESULT} 0 PARENT_SCOPE )
    ENDIF()
ENDFUNCTION()


# Add a provisional test
FUNCTION( ADD_AMP_PROVISIONAL_TEST EXEFILE )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Check if test has already been added
    GET_TARGET_PROPERTY(tmp ${EXEFILE} LOCATION)
    IF ( NOT tmp )
        # The target has not been added
        SET( CXXFILE ${EXEFILE}.cc )
        SET( TESTS_SO_FAR ${TESTS_SO_FAR} ${EXEFILE} )
        IF ( NOT EXCLUDE_TESTS_FROM_ALL )
            ADD_EXECUTABLE( ${EXEFILE} ${CXXFILE} )
        ELSE()
            ADD_EXECUTABLE( ${EXEFILE} EXCLUDE_FROM_ALL ${CXXFILE} )
        ENDIF()
        ADD_AMP_EXE_DEP( ${EXEFILE} )
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(Configuration)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSEIF( ${tmp} STREQUAL "${CMAKE_CURRENT_BINARY_DIR}/$(OutDir)/${EXEFILE}.exe" )
        # The correct target has already been added
    ELSE()
        # We are trying to add 2 different tests with the same name
        MESSAGE ( "Existing test: ${tmp}" )
        MESSAGE ( "New test:      ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE}" )
        MESSAGE ( FATAL_ERROR "Trying to add 2 different tests with the same name" )
    ENDIF()
ENDFUNCTION()


# Macro to create the test name
MACRO( CREATE_TEST_NAME TEST ${ARGN} )
    IF ( PACKAGE )
        SET( TESTNAME "${PACKAGE}::${TEST}" )
    ELSE()
        SET( TESTNAME "${TEST}" )
    ENDIF()
    foreach( tmp ${ARGN})
        SET( TESTNAME "${TESTNAME}--${tmp}")
    endforeach()
    # STRING(REGEX REPLACE "--" "-" TESTNAME ${TESTNAME} )
ENDMACRO()


# Add a executable as a test
FUNCTION( ADD_AMP_TEST EXEFILE ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_AMP_PROVISIONAL_TEST ( ${EXEFILE} )
    CREATE_TEST_NAME( ${EXEFILE} ${ARGN} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_EXT_MPI_FOR_SERIAL_TESTS )
        ADD_TEST ( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
    ELSE()
        ADD_TEST ( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS 1 )
ENDFUNCTION()


# Add a executable as a weekly test
FUNCTION( ADD_AMP_WEEKLY_TEST EXEFILE PROCS ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_AMP_PROVISIONAL_TEST ( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF( ${PROCS} STREQUAL "1" )
        CREATE_TEST_NAME( "${EXEFILE}_WEEKLY" ${ARGN} )
        IF( USE_EXT_MPI_FOR_SERIAL_TESTS )
            ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${EXE} ${ARGN} )
        ELSE()
            ADD_TEST( ${TESTNAME} ${CMAKE_CURRENT_BINARY_DIR}/${EXEFILE} ${ARGN} )
        ENDIF()
    ELSEIF( USE_EXT_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs_WEEKLY" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
    ENDIF()
    SET_TESTS_PROPERTIES ( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${PROCS} )
ENDFUNCTION()


# Add a executable as a parallel test
FUNCTION( ADD_AMP_TEST_PARALLEL EXEFILE PROCS ${ARGN} )
    # Check if we actually want to add the test
    KEEP_TEST( RESULT )
    IF ( NOT RESULT )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_AMP_PROVISIONAL_TEST( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    IF ( USE_EXT_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "${EXEFILE}_${PROCS}procs" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${PROCS} )
    ENDIF()
ENDFUNCTION()


# Add a executable as an example
FUNCTION( ADD_AMP_EXAMPLE EXEFILE PROCS ${ARGN} )
    # Add the file to the example doxygen file
    SET( VALUE 0 )
    FOREACH(_variableName ${EXAMPLE_LIST})
        IF ( "${_variableName}" STREQUAL "${EXEFILE}" )
            SET( VALUE 1 )
        ENDIF()
    ENDFOREACH()
    IF ( NOT ${VALUE} )
        FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "* \\ref ${EXEFILE} \"${EXEFILE}\"\n" )
        SET( EXAMPLE_LIST ${EXAMPLE_LIST} ${EXEFILE} CACHE INTERNAL "example_list" FORCE )
    ENDIF()
    # Check if we actually want to add the test
    IF ( ONLY_BUILD_DOCS )
        RETURN()
    ENDIF()
    # Add the provisional test
    ADD_AMP_PROVISIONAL_TEST( ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE ${EXEFILE} LOCATION)
    STRING(REGEX REPLACE "\\$\\(Configuration\\)" "${CONFIGURATION}" EXE "${EXE}" )
    ADD_DEPENDENCIES( build-examples ${EXEFILE} )
    GET_TARGET_PROPERTY(EXE2 ${EXEFILE} LOCATION)
    ADD_CUSTOM_COMMAND( TARGET ${EXEFILE} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${EXE2} "${EXAMPLE_INSTALL_DIR}/${EXEFILE}"
    )
    IF( ${PROCS} STREQUAL "1" AND (NOT USE_EXT_MPI_FOR_SERIAL_TESTS) )
        CREATE_TEST_NAME( "example--${EXEFILE}" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${EXE} ${ARGN} )    
    ELSEIF ( USE_EXT_MPI AND NOT (${PROCS} GREATER ${TEST_MAX_PROCS}) )
        CREATE_TEST_NAME( "example--${EXEFILE}_${PROCS}procs" ${ARGN} )
        ADD_TEST( ${TESTNAME} ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${EXE} ${ARGN} )
        SET_TESTS_PROPERTIES( ${TESTNAME} PROPERTIES FAIL_REGULAR_EXPRESSION ".*FAILED.*" PROCESSORS ${PROCS} )
    ENDIF()
ENDFUNCTION()


# Begin configure for the examples for a package
MACRO( BEGIN_EXAMPLE_CONFIG PACKAGE )
    # Set example install dir
    SET( EXAMPLE_INSTALL_DIR ${AMP_INSTALL_DIR}/examples/${PACKAGE} )
    # Create list of examples
    SET( EXAMPLE_LIST "dummy" CACHE INTERNAL "example_list" FORCE )
    # Create doxygen input file for examples
    SET( DOXYFILE_EXTRA_SOURCES ${DOXYFILE_EXTRA_SOURCES} ${EXAMPLE_INSTALL_DIR} CACHE INTERNAL "doxyfile_extra_sources") 
    FILE(WRITE  ${EXAMPLE_INSTALL_DIR}/examples.h "// Include file for doxygen providing the examples for ${PACKAGE}\n")
    FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "/*! \\page Examples_${PACKAGE}\n" )
ENDMACRO()

# Install the examples
MACRO( INSTALL_AMP_EXAMPLE PACKAGE )
    FILE(APPEND ${EXAMPLE_INSTALL_DIR}/examples.h "*/\n" )
    SET( EXAMPLE_INSTALL_DIR "" )
ENDMACRO()


# Add a test with 1, 2, and 4 processors
MACRO( ADD_AMP_TEST_1_2_4 EXENAME ${ARGN} )
    ADD_AMP_TEST ( ${EXENAME} ${ARGN} )
    ADD_AMP_TEST_PARALLEL ( ${EXENAME} 2 ${ARGN} )
    ADD_AMP_TEST_PARALLEL ( ${EXENAME} 4 ${ARGN} )
ENDMACRO()


# Macro to check if a flag is enabled
MACRO ( CHECK_ENABLE_FLAG FLAG DEFAULT )
    IF( NOT DEFINED ${FLAG} )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ${FLAG}  STREQUAL "" )
        SET( ${FLAG} ${DEFAULT} )
    ELSEIF( ( ${${FLAG}} STREQUAL "false" ) OR ( ${${FLAG}} STREQUAL "0" ) OR ( ${${FLAG}} STREQUAL "OFF" ) )
        SET( ${FLAG} 0 )
    ELSEIF( ( ${${FLAG}} STREQUAL "true" ) OR ( ${${FLAG}} STREQUAL "1" ) OR ( ${${FLAG}} STREQUAL "ON" ) )
        SET( ${FLAG} 1 )
    ELSE()
        MESSAGE( "Bad value for ${FLAG} (${${FLAG}}); use true or false" )
    ENDIF ()
ENDMACRO()


# Macro to check if a compiler flag is valid
MACRO (CHECK_C_COMPILER_FLAG _FLAG _RESULT)
    SET(SAFE_CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS}")
    SET(CMAKE_REQUIRED_DEFINITIONS "${_FLAG}")
    CHECK_C_SOURCE_COMPILES("int main() { return 0;}" ${_RESULT}
        # Some compilers do not fail with a bad flag
        FAIL_REGEX "error: bad value (.*) for .* switch"       # GNU
        FAIL_REGEX "argument unused during compilation"        # clang
        FAIL_REGEX "is valid for .* but not for C"             # GNU
        FAIL_REGEX "unrecognized .*option"                     # GNU
        FAIL_REGEX "ignoring unknown option"                   # MSVC
        FAIL_REGEX "[Uu]nknown option"                         # HP
        FAIL_REGEX "[Ww]arning: [Oo]ption"                     # SunPro
        FAIL_REGEX "command option .* is not recognized"       # XL
        FAIL_REGEX "WARNING: unknown flag:"                    # Open64
        FAIL_REGEX " #10159: "                                 # ICC
    )
    SET(CMAKE_REQUIRED_DEFINITIONS "${SAFE_CMAKE_REQUIRED_DEFINITIONS}")
ENDMACRO(CHECK_C_COMPILER_FLAG)


# Macro to print all variables
MACRO( PRINT_ALL_VARIABLES )
    GET_CMAKE_PROPERTY(_variableNames VARIABLES)
    FOREACH(_variableName ${_variableNames})
        message(STATUS "${_variableName}=${${_variableName}}")
    ENDFOREACH()
ENDMACRO()


# Macro to change the classification of a package
MACRO( SET_PACKAGE_CLASSIFICATION  PACKAGE_LIST  PACKAGE_NAME  CLASS )
    LIST(FIND ${PACKAGE_LIST} ${PACKAGE_NAME} PACKAGE_NAME_IDX)
    IF (PACKAGE_NAME_IDX EQUAL -1)
        MESSAGE(FATAL_ERROR "Package ${PACKAGE_NAME} not found in list of packages!")
    ELSE()
        MATH(EXPR PACKAGE_CLASSIFICATION_IDX "${PACKAGE_NAME_IDX}+2")
        LIST(INSERT ${PACKAGE_LIST} ${PACKAGE_CLASSIFICATION_IDX} ${CLASS})
        MATH(EXPR PACKAGE_CLASSIFICATION_IDX "${PACKAGE_CLASSIFICATION_IDX} + 1")
        LIST (REMOVE_AT ${PACKAGE_LIST} ${PACKAGE_CLASSIFICATION_IDX})
    ENDIF()
ENDMACRO()


# Macro to "disable" a package on the given platform (this mearly changes it to experimental)
MACRO( PACKAGE_DISABLE_ON_PLATFORMS  PACKAGE_LIST  PACKAGE_NAME )
    FOREACH(HOSTTYPE ${ARGN})
        IF (${PROJECT_NAME}_HOSTTYPE STREQUAL ${HOSTTYPE})
            SET_PACKAGE_CLASSIFICATION(${PACKAGE_LIST} ${PACKAGE_NAME} EX)
            IF (${PROJECT_NAME}_ENABLE_${PACKAGE_NAME})
                MESSAGE(
                  "\n***"
                  "\n*** WARNING: User has set ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}=ON but the"
                  "\n*** package ${PACKAGE_NAME} is not supported on this platform type '${HOSTTYPE}'!"
                  "\n***\n"
               )
            ENDIF()
        ENDIF()
    ENDFOREACH()
ENDMACRO()


# add custom target distclean
# cleans and removes cmake generated files etc.
MACRO ( ADD_DISTCLEAN )
    SET(DISTCLEANED
        cmake.depends
        cmake.check_depends
        CMakeCache.txt
        CMakeFiles
        CMakeTmp
        cmake.check_cache
        *.cmake
        compile.log
        Doxyfile
        Makefile
        core core.*
        src
        ampdir
        AMP
        DartConfiguration.tcl
        Testing
        install_manifest.txt
        nek
    )
    ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
    IF (UNIX)
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND rm
            ARGS    -Rf ${DISTCLEANED}
            TARGET  distclean
        )
    ELSE()
        SET( DISTCLEANED
            ${DISTCLEANED}
            *.vcxproj*
            ipch
            x64
        )
        FILE(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "del /s /q /f " )
        FOREACH (fileToDelete ${DISTCLEANED})
            FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "${fileToDelete} " )
        ENDFOREACH ()
        FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "\n" )
        FOREACH (fileToDelete ${DISTCLEANED})
            FILE(APPEND ${CMAKE_CURRENT_BINARY_DIR}/distclean.bat "for /d %%x in (${fileToDelete}) do rd /s /q \"%%x\"\n" )
        ENDFOREACH ()
        ADD_CUSTOM_COMMAND(
            DEPENDS clean
            COMMENT "distribution clean"
            COMMAND distclean.bat & del /s/q/f distclean.bat
            TARGET  distclean
        )
    ENDIF()
ENDMACRO()


# Save the necessary cmake variables to a file for applications to load
# Note: we need to save the external libraries in the same order as AMP for consistency
MACRO ( SAVE_CMAKE_FLAGS )
    # Write the header (comments)
    file(WRITE  ${AMP_INSTALL_DIR}/amp.cmake "# This is a automatically generate file to include AMP within another application\n\n" )
    # Write the compilers and compile flags
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Set the compilers and compile flags\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_C_COMPILER ${CMAKE_C_COMPILER})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_Fortran_COMPILER ${CMAKE_Fortran_COMPILER})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_FORTRAN ${USE_FORTRAN})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_C_FLAGS \"${CMAKE_C_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(LDFLAGS \"${LDFLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(COMPILE_MODE ${COMPILE_MODE})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS \"${CMAKE_EXE_LINK_DYNAMIC_C_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS \"${CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_SHARED_LIBRARY_C_FLAGS \"${CMAKE_SHARED_LIBRARY_C_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS \"${CMAKE_SHARED_LIBRARY_CXX_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_SHARED_LINKER_FLAGS \"${CMAKE_SHARED_LINKER_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS \"${CMAKE_SHARED_LIBRARY_LINK_C_FLAGS}\")\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS \"${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS}\")\n" )
    # Write the AMP_DATA and AMP_SOURCE paths
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Set the AMP data and source directories\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(AMP_DATA ${AMP_DATA})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(AMP_SOURCE ${AMP_SOURCE_DIR})\n" )
    # Create the AMP libraries and include paths
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Set the AMP libraries\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(AMP_LIBS ${AMP_LIBS})\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES ( ${AMP_TRUNK}/external/boost/include )\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES ( ${AMP_INSTALL_DIR}/include )\n" )
    IF ( USE_AMP_UTILS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_UTILS ${USE_AMP_UTILS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_UTILS ) \n" )
    ENDIF()
    IF ( USE_AMP_MESH )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_MESH ${USE_AMP_MESH}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_MESH ) \n" )
    ENDIF()
    IF ( USE_AMP_DISCRETIZATION )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_DISCRETIZATION ${USE_AMP_DISCRETIZATION}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_DISCRETIZATION ) \n" )
    ENDIF()
    IF ( USE_AMP_VECTORS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_VECTORS ${USE_AMP_VECTORS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_VECTORS ) \n" )
    ENDIF()
    IF ( USE_AMP_MATRICES )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_MATRICES ${USE_AMP_MATRICES}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_MATRICES ) \n" )
    ENDIF()
    IF ( USE_AMP_MATERIALS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_MATERIALS ${USE_AMP_MATERIALS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_MATERIALS ) \n" )
    ENDIF()
    IF ( USE_AMP_OPERATORS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_OPERATORS ${USE_AMP_OPERATORS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_OPERATORS ) \n" )
    ENDIF()
    IF ( USE_AMP_SOLVERS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_SOLVERS ${USE_AMP_SOLVERS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_SOLVERS ) \n" )
    ENDIF()
    IF ( USE_AMP_TIME_INTEGRATORS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_AMP_TIME_INTEGRATORS ${USE_AMP_TIME_INTEGRATORS}) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_AMP_TIME_INTEGRATORS ) \n" )
    ENDIF()
    # Create the external libraries and include paths in the order they are linked in AMP
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS )\n" )
    # Add boost
    IF ( USE_EXT_BOOST )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add boost\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${BOOST_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_BOOST 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_BOOST ) \n" )
    ENDIF()
    # Add Libmesh
    IF ( USE_EXT_LIBMESH )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add Libmesh\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${LIBMESH_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${LIBMESH_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -DLIBMESH_ENABLE_PARMESH )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_LIBMESH 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_LIBMESH ) \n" )
    ENDIF()
    # Add NEK
    IF ( USE_EXT_NEK )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add NEK\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${NEK_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${NEK_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_NEK 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_NEK ) \n" )
    ENDIF()
    # Add MOAB
    IF ( USE_EXT_MOAB )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add MOAB\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${MOAB_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${MOAB_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_MOAB 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_MOAB ) \n" )
    ENDIF()
    # Add DENDRO
    IF ( USE_EXT_DENDRO )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add DENDRO\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${DENDRO_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${DENDRO_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_DENDRO 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_DENDRO ) \n" )
    ENDIF()
    # Add Netcdf
    IF ( USE_NETCDF )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add Netcdf\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${NETCDF_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${NETCDF_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_NETCDF 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_NETCDF ) \n" )
    ENDIF()
    # Add Trilinos
    IF ( USE_EXT_TRILINOS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add Trilinos\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${TRILINOS_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${TRILINOS_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_TRILINOS 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_UTILS   ${USE_TRILINOS_UTILS}   ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_TEUCHOS ${USE_TRILINOS_TEUCHOS} ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_VECTORS ${USE_TRILINOS_VECTORS} ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_SOLVERS ${USE_TRILINOS_SOLVERS} ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_THYRA   ${USE_TRILINOS_THYRA}   ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_NOX     ${USE_TRILINOS_NOX}     ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_STKMESH ${USE_TRILINOS_STKMESH} ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_TRILINOS_STRATIMIKOS ${USE_TRILINOS_STRATIMIKOS} ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_TRILINOS ) \n" )
        IF ( USE_TRILINOS_THYRA )
            file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_TRILINOS_THYRA ) \n" )
        ENDIF()
        IF ( USE_TRILINOS_NOX )
            file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_TRILINOS_NOX ) \n" )
        ENDIF()
        IF ( USE_TRILINOS_STKMESH )
            file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_TRILINOS_STKMESH ) \n" )
        ENDIF()
    ENDIF()
    # Add PETsc
    IF ( USE_EXT_PETSC )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add PETsc\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${PETSC_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${PETSC_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_PETSC 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_PETSC ) \n" )
    ENDIF()
    # Add Sundials
    IF ( USE_EXT_SUNDIALS )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add Sundials\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${SUNDIALS_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${SUNDIALS_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_SUNDIALS 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_SUNDIALS ) \n" )
    ENDIF()
    # Add Silo
    IF ( USE_EXT_SILO )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add silo\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${SILO_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${SILO_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_SILO 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_SILO ) \n" )
    ENDIF()
    # Add Hypre
    IF ( USE_EXT_HYPRE )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add hypre\n" )
        # file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${HYPRE_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${HYPRE_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_HYPRE 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_HYPRE ) \n" )
    ENDIF()
    # Add X11
    IF ( USE_EXT_X11 )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add X11\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${X11_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${X11_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_X11 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_X11 ) \n" )
    ENDIF()
    # Add HDF5
    IF ( USE_EXT_HDF5 )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add HDF5\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${HDF5_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${HDF5_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_HDF5 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_HDF5 ) \n" )
    ENDIF()
    # Add MPI
    IF ( USE_EXT_MPI )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add MPI\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "INCLUDE_DIRECTORIES( ${MPI_INCLUDE} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS}  ${MPI_LINK_FLAGS} ${MPI_LIBRARIES} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_MPI 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -D USE_EXT_MPI ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(MPIEXEC ${MPIEXEC} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_EXT_MPI_FOR_SERIAL_TESTS ${USE_EXT_MPI_FOR_SERIAL_TESTS} )\n" )
    ELSE()
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET(USE_EXT_MPI_FOR_SERIAL_TESTS 0 )\n" )
    ENDIF()
    # Add LAPACK and BLAS
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add LAPACK/BLAS\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS}  ${LAPACK_LIBS} ${BLAS_LIBS} ${BLAS_LAPACK_LIBS} )\n" )
    # Add LDLIBS
    IF ( LDLIBS )
    	file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS}  ${LDLIBS} )\n" )
    ENDIF()
    # Add coverage
    IF ( ENABLE_GCOV )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add coverage flags\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${COVERAGE_LIBS} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( ENABLE_GCOV 1 ) \n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "ADD_DEFINITIONS( -fprofile-arcs -ftest-coverage ) \n" )
    ENDIF ()
    # Add doxygen
    IF ( USE_EXT_DOXYGEN )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add doxygen flags\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( USE_EXT_DOXYGEN ${USE_EXT_DOXYGEN} )\n" )
        file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( DOXYGEN_MACROS \"${DOXYGEN_MACROS}\" )\n" )
    ENDIF ()
    # Add misc flags
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "# Add misc flags\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( EXTERNAL_LIBS $" "{EXTERNAL_LIBS} ${SYSTEM_LIBS} )\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "SET( TEST_MAX_PROCS ${TEST_MAX_PROCS} )\n" )
    file(APPEND ${AMP_INSTALL_DIR}/amp.cmake "\n" )
ENDMACRO ()


# Add an external subdirectory
MACRO ( ADD_EXTERNAL_PACKAGE_SUBDIRECTORY SUBDIR_NAME SUBDIR_PATH )
    VERIFY_PATH ( ${SUBDIR_PATH} )
    FIND_FILES_PATH ( ${SUBDIR_PATH} )
    FILE ( GLOB HFILES RELATIVE ${SUBDIR_PATH} ${SUBDIR_PATH}/*.h ${SUBDIR_PATH}/*.hh ${SUBDIR_PATH}/*.I )
    FOREACH (HFILE ${HFILES})
        CONFIGURE_FILE ( ${SUBDIR_PATH}/${HFILE} ${AMP_INSTALL_DIR}/include/${CURPACKAGE}/${HFILE} COPYONLY )
    ENDFOREACH ()
    ADD_SUBDIRECTORY ( ${SUBDIR_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${SUBDIR_NAME} )
ENDMACRO ()


