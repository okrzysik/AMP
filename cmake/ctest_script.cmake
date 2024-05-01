# ctest script for building, running, and submitting the test results 
# Usage:  ctest -S script,build
#   build = debug / optimized / valgrind
# Note: this test will use use the number of processors defined in the variable N_PROCS,
#   the enviornmental variable N_PROCS, or the number of processors availible (if not specified)


# Set the Project variables
SET( PROJ AMP )


# Set platform specific variables
SITE_NAME( HOSTNAME )
STRING( REGEX REPLACE  "-(ext|login)(..|.)"  ""  HOSTNAME  "${HOSTNAME}" )
SET( MPIEXEC             $ENV{MPIEXEC}            )
SET( TPL_DIRECTORY      "$ENV{TPL_DIRECTORY}"     )
SET( AMP_DATA            $ENV{AMP_DATA}           )
SET( DISABLE_GXX_DEBUG   $ENV{DISABLE_GXX_DEBUG}  )
SET( COVERAGE_COMMAND    $ENV{COVERAGE_COMMAND}   )
SET( VALGRIND_COMMAND    $ENV{VALGRIND_COMMAND}   )
SET( CMAKE_MAKE_PROGRAM  $ENV{CMAKE_MAKE_PROGRAM} )
SET( CTEST_CMAKE_GENERATOR $ENV{CTEST_CMAKE_GENERATOR} )
SET( BUILD_SERIAL        $ENV{BUILD_SERIAL}       )
SET( DISABLE_FORTRAN     $ENV{DISABLE_FORTRAN}    )
SET( SKIP_TESTS          $ENV{SKIP_TESTS}         )
SET( BUILDNAME_POSTFIX  "$ENV{BUILDNAME_POSTFIX}" )
SET( TPL_LIST_REQUIRED  "$ENV{TPL_LIST_REQUIRED}" )
SET( TPL_LIST_OPTIONAL  "$ENV{TPL_LIST_OPTIONAL}" )
SET( USE_CUDA            $ENV{USE_CUDA}           )
SET( CUDA_FLAGS         "$ENV{CUDA_FLAGS}"        )
SET( CTEST_SITE         "$ENV{CTEST_SITE}"        )
SET( CTEST_URL          "$ENV{CTEST_URL}"         )
SET( EXCLUDE_WEEKLY     "$ENV{EXCLUDE_WEEKLY}"    )


# Get the source directory based on the current directory
IF ( NOT ${PROJ}_SOURCE_DIR )
    SET( ${PROJ}_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." )
ENDIF()
IF ( NOT CMAKE_MAKE_PROGRAM )
    SET( CMAKE_MAKE_PROGRAM make )
ENDIF()


# Check that we specified the build type to run
IF( NOT CTEST_SCRIPT_ARG )
    MESSAGE(FATAL_ERROR "No build specified: ctest -S /path/to/script,build (debug/optimized/valgrind")
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "debug" )
    SET( CTEST_BUILD_NAME "${PROJ}-debug" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND ${COVERAGE_COMMAND} )
    SET( ENABLE_GCOV "true" )
    SET( USE_VALGRIND FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "optimized") OR (${CTEST_SCRIPT_ARG} STREQUAL "opt") )
    SET( CTEST_BUILD_NAME "${PROJ}-opt" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
ELSEIF( (${CTEST_SCRIPT_ARG} STREQUAL "weekly") )
    SET( CTEST_BUILD_NAME "AMP-Weekly" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND FALSE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "valgrind" )
    SET( CTEST_BUILD_NAME "${PROJ}-valgrind" )
    SET( CMAKE_BUILD_TYPE "Debug" )
    SET( CTEST_COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
    SET( USE_VALGRIND TRUE )
ELSEIF( ${CTEST_SCRIPT_ARG} STREQUAL "doc" )
    SET( CTEST_BUILD_NAME "${PROJ}-doc" )
    SET( CMAKE_BUILD_TYPE "Release" )
    SET( BUILD_ONLY_DOCS "true" )
ELSE()
    MESSAGE(FATAL_ERROR "Invalid build (${CTEST_SCRIPT_ARG}): ctest -S /path/to/script,build (debug/opt/valgrind")
ENDIF()
IF ( BUILDNAME_POSTFIX )
    SET( CTEST_BUILD_NAME "${CTEST_BUILD_NAME}-${BUILDNAME_POSTFIX}" )
ENDIF()
IF ( NOT COVERAGE_COMMAND )
    SET( ENABLE_GCOV "false" )
ENDIF()


# Set the number of processors
SET( N_PROCS $ENV{N_PROCS} )
IF ( NOT DEFINED N_PROCS )
    SET( N_PROCS 4 ) # Default number of processor if all else fails
    IF ( EXISTS "/proc/cpuinfo" )
        # Linux
        FILE( STRINGS "/proc/cpuinfo" procs REGEX "^processor.: [0-9]+$" )
        LIST( LENGTH procs N_PROCS )
    ELSEIF( APPLE )
        FIND_PROGRAM( cmd_sys_pro "system_profiler" )
        IF ( cmd_sys_pro )
            EXECUTE_PROCESS( COMMAND ${cmd_sys_pro} OUTPUT_VARIABLE info )
            STRING( REGEX REPLACE "^.*Total Number of Cores: ([0-9]+).*$" "\\1" N_PROCS "${info}" )
        ENDIF()
    ENDIF()
ENDIF()


# Set the nightly start time
# This controls the version of a checkout from cvs/svn (ignored for mecurial/git)
# This does not control the start of the day displayed on CDash, that is controled by the CDash project settings
SET( NIGHTLY_START_TIME "$ENV{NIGHTLY_START_TIME}" )
IF ( NOT NIGHTLY_START_TIME )
    SET( NIGHTLY_START_TIME "17:00:00 EST" )
ENDIF()
SET( CTEST_NIGHTLY_START_TIME ${NIGHTLY_START_TIME} )


# Set basic variables
SET( CTEST_PROJECT_NAME "${PROJ}" )
SET( CTEST_SOURCE_DIRECTORY "${${PROJ}_SOURCE_DIR}" )
SET( CTEST_BINARY_DIRECTORY "." )
SET( CTEST_DASHBOARD "Nightly" )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500 )
SET( CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500 )
SET( CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE 10000 )
SET( CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE 20000 )
SET( CTEST_COMMAND "\"${CTEST_EXECUTABLE_NAME}\" -D ${CTEST_DASHBOARD}" )
SET( CTEST_BUILD_COMMAND "${CMAKE_MAKE_PROGRAM} -i -j ${N_PROCS} install" )
SET( CTEST_CUSTOM_WARNING_EXCEPTION 
    "has no symbols"
    "the table of contents is empty"
    "warning: -jN forced in submake: disabling jobserver mode" 
    "warning: jobserver unavailable" 
    "This object file does not define any previously undefined public symbols"
)


# Set timeouts: 10 minutes for debug, 5 for opt, and 30 minutes for valgrind/weekly
IF ( USE_VALGRIND )
    SET( CTEST_TEST_TIMEOUT 1800 )
ELSEIF( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    SET( CTEST_TEST_TIMEOUT 600 )
ELSE()
    SET( CTEST_TEST_TIMEOUT 300 )
ENDIF()


# Set valgrind options
#SET (VALGRIND_COMMAND_OPTIONS "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --track-origins=yes --malloc-fill=0xff --free-fill=0xfe --suppressions=${AMP_SOURCE_DIR}/ValgrindSuppresionFile" )
SET( VALGRIND_COMMAND_OPTIONS  "--tool=memcheck --leak-check=yes --track-fds=yes --num-callers=50 --show-reachable=yes --suppressions=${AMP_SOURCE_DIR}/ValgrindSuppresionFile" )
IF ( USE_VALGRIND )
    SET( MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECKCOMMAND ${VALGRIND_COMMAND} )
    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS ${VALGRIND_COMMAND_OPTIONS} )
    SET( CTEST_MEMORYCHECKCOMMAND_OPTIONS  ${VALGRIND_COMMAND_OPTIONS} )
ENDIF()


# Clear the binary directory and create an initial cache
EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt )
EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles )
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "CTEST_TEST_CTEST:BOOL=1")


# Set the configure options
SET( CTEST_OPTIONS "-DAMP_DATA='${AMP_DATA}'" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DTPL_DIRECTORY='${TPL_DIRECTORY}'" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DAMP_INSTALL_DIR='${CMAKE_CURRENT_BINARY_DIR}/ampdir'" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DEXCLUDE_TESTS_FROM_ALL:BOOL=FALSE" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DLDLIBS:STRING='${LDLIBS}'" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DENABLE_GCOV:BOOL=${ENABLE_GCOV}" )
SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_EXT_DOXYGEN=0" )
IF ( USE_CUDA )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_CUDA=1;-DCMAKE_CUDA_FLAGS=${CUDA_FLAGS}" )
ELSE()
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_CUDA=0" )
ENDIF()
IF ( DISABLE_FORTRAN )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DUSE_FORTRAN=OFF;-DUSE_EXT_FORTRAN=OFF" )
ENDIF()
IF ( TPL_LIST_REQUIRED )
    STRING( REGEX REPLACE  ";"  ","  TPL_LIST_REQUIRED  "${TPL_LIST_REQUIRED}" )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DTPL_LIST_REQUIRED='${TPL_LIST_REQUIRED}'" )
ENDIF()
IF ( TPL_LIST_OPTIONAL )
    STRING( REGEX REPLACE  ";"  ","  TPL_LIST_OPTIONAL  "${TPL_LIST_OPTIONAL}" )
    SET( CTEST_OPTIONS "${CTEST_OPTIONS};-DTPL_LIST_OPTIONAL='${TPL_LIST_OPTIONAL}'" )
ENDIF()
MESSAGE("Configure options:")
MESSAGE("   ${CTEST_OPTIONS}")


# Configure the drop site
IF ( NOT CTEST_SITE )
    SET( CTEST_SITE ${HOSTNAME} )
ENDIF()
IF ( NOT CTEST_URL )
    SET( CTEST_DROP_METHOD "http" )
    SET( CTEST_DROP_LOCATION "/CDash/submit.php?project=AMP" )
    SET( CTEST_DROP_SITE_CDASH TRUE )
    SET( DROP_SITE_CDASH TRUE )
    SET( CTEST_DROP_SITE ${CTEST_SITE} )
ELSE()
    STRING( REPLACE "PROJECT" "AMP" CTEST_URL "${CTEST_URL}" )
    SET( CTEST_SUBMIT_URL "${CTEST_URL}" )
ENDIF()


# Configure and run the tests
CTEST_START( "${CTEST_DASHBOARD}" )
CTEST_UPDATE()
CTEST_SUBMIT( PARTS Update )
CTEST_CONFIGURE(
    BUILD   ${CTEST_BINARY_DIRECTORY}
    SOURCE  ${CTEST_SOURCE_DIRECTORY}
    OPTIONS "${CTEST_OPTIONS}"
)
CTEST_SUBMIT( PARTS Configure )


# Run the configure/build/test
CTEST_BUILD()
CTEST_SUBMIT( PARTS Build )
EXECUTE_PROCESS( COMMAND ${CMAKE_MAKE_PROGRAM} install )
IF ( SKIP_TESTS )
    # Do not run tests
    SET( CTEST_COVERAGE_COMMAND )
ELSEIF ( USE_VALGRIND )
    CTEST_MEMCHECK( EXCLUDE "(procs|WEEKLY|cppcheck|cppclean|test_crash)"  PARALLEL_LEVEL ${N_PROCS} )
ELSEIF ( EXCLUDE_WEEKLY )
    CTEST_TEST( EXCLUDE WEEKLY  PARALLEL_LEVEL ${N_PROCS} )
ELSE()
    CTEST_TEST( PARALLEL_LEVEL ${N_PROCS} )
ENDIF()
IF( CTEST_COVERAGE_COMMAND )
    CTEST_COVERAGE()
ENDIF()
CTEST_SUBMIT( PARTS Test )
CTEST_SUBMIT( PARTS Coverage )
CTEST_SUBMIT( PARTS MemCheck )
CTEST_SUBMIT( PARTS Done )


# Write a message to test for success in the ctest-builder
MESSAGE( "ctest_script ran to completion" )


