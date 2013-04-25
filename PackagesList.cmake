INCLUDE(TribitsListHelpers)
#
# Define the package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every push to the global repo.  The
#     default enable for PS packages is empty "" which allows the PS package
#     to be enabled implicitly based on other criteria.  The option
#     AMP_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every push to the
#     global repo.  For example, a package that does not provide any
#     significant functionally is classified as a SS package even
#     through it could be classified as PS just based on its required package
#     and TPL dependencies.  SS packages will have their default enables set
#     to empty "".  This allows them to be enabled implicilty.  When
#     AMP_ENABLE_ALL_PACKAGES=ON but AMP_ENABLE_SECONDARY_STABLE_CODE=OFF, 
#     the SS packages will not be enabled.  However, when AMP_ENABLE_ALL_PACKAGES=ON 
#     and AMP_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in pre-push testing.  However, SS packages are
#     tested by the post-push CI and nightly testing processes.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in pre-push testring and are not tested as part of the
#     post-push CI or nightly testing processes.  However, package developers
#     of EX pacakges are encouraged to set up their own nightly testing for
#     thier EX packages.
#
# NOTE: These packages must be listed in strictly ascending order in terms of
# package dependencies.  If you get the order wrong, then an error message
# will be printed during configuration with CMake.
#

CMAKE_POLICY(SET CMP0014 OLD)

SET( AMP_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    AMP_UTILITIES           src/utils               PS  # Does not have any required dependencies
    AMP_MESH                src/ampmesh             PS  # Requires utils
    AMP_DISCRETIZATION      src/discretization      PS  # Requires utils, ampmesh
    AMP_VECTORS             src/vectors             PS  # Requires utils
    AMP_MATRICES            src/matrices            PS  # Requires utils
    AMP_MATERIALS           src/materials           PS  # Requires utils
    AMP_OPERATORS           src/operators           SS  # Requires utils, ampmesh, discretization, vectors, matrices, materials, libmesh
    AMP_TIME_INTEGRATORS    src/time_integrators    SS  # Requires utils, ampmesh, discretization, vectors, matrices, materials, operators
    AMP_SOLVERS             src/solvers             SS  # Requires utils, ampmesh, discretization, vectors, matrices, materials, operators
)


#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#
PACKAGE_DISABLE_ON_PLATFORMS(operators Windows)
PACKAGE_DISABLE_ON_PLATFORMS(time_integrators Windows)
PACKAGE_DISABLE_ON_PLATFORMS(solvers Windows)


# Set the AMP source, build, and install directories
IF ( NOT AMP_SOURCE_DIR )
    STRING(REGEX REPLACE "/PackagesList.cmake" "" AMP_SOURCE_DIR ${CMAKE_CURRENT_LIST_FILE} )
ENDIF()
SET (AMP_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} )
IF ( NOT AMP_INSTALL_DIR )
  SET ( AMP_INSTALL_DIR ${AMP_BUILD_DIR}/ampdir )
ENDIF()
MESSAGE ( "Installing AMP in "${AMP_INSTALL_DIR} )

# Initialize the libaries (flags will be overwritten when the libraries are configured)
SET( USE_EXT_MPI 0 )
SET( TEST_MAX_PROCS 8 )

# Create custom targets for build-test, check, and distclean
INCLUDE ( ${AMP_SOURCE_DIR}/cmake/macros.cmake )
INCLUDE ( ${AMP_SOURCE_DIR}/cmake/libraries.cmake )
ADD_CUSTOM_TARGET ( build-test )
ADD_CUSTOM_TARGET ( check COMMAND  make test  )
ADD_DISTCLEAN()

# Check if we want to enable fortran
ENABLE_LANGUAGE(C)
ENABLE_LANGUAGE(CXX)
IF ( DEFINED USE_EXT_FORTRAN )
    SET ( USE_FORTRAN ${USE_EXT_FORTRAN} )
ENDIF()
IF ( NOT DEFINED USE_FORTRAN )
    SET ( USE_FORTRAN 1 )
ELSEIF ( ( ${USE_FORTRAN} STREQUAL "false" ) OR ( ${USE_FORTRAN} STREQUAL "0" ) OR ( ${USE_FORTRAN} STREQUAL "OFF" ) )
    SET ( USE_FORTRAN 0 )
ELSEIF ( ( ${USE_FORTRAN} STREQUAL "true" ) OR ( ${USE_FORTRAN} STREQUAL "1" ) OR ( ${USE_FORTRAN} STREQUAL "ON" ) )
    SET (USE_FORTRAN 1 )
ELSE()
    MESSAGE ( "Bad value for USE_FORTRAN; use true or false" )
ENDIF()
IF ( USE_FORTRAN )
    ENABLE_LANGUAGE (Fortran)
ENDIF()

# Set the compile flags
CONFIGURE_SYSTEM()
IF ( CMAKE_BUILD_TYPE ) 
    IF ( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR ${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
        SET_DEBUG_MACROS()
    ELSEIF ( ${CMAKE_BUILD_TYPE} STREQUAL "Release" OR ${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
        SET_OPTIMIZED_MACROS()
    ELSE()
        MESSAGE(FATAL_ERROR "Unkown value for CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
    ENDIF()
ELSE()
    VERIFY_VARIABLE ( "COMPILE_MODE" )
    IF ( ${COMPILE_MODE} STREQUAL "debug" )
        SET(CMAKE_BUILD_TYPE "Debug")
        SET_DEBUG_MACROS()
    ELSEIF ( ${COMPILE_MODE} STREQUAL "optimized" )
        SET(CMAKE_BUILD_TYPE "Release")
        SET_OPTIMIZED_MACROS()
    ELSE()
        MESSAGE ( FATAL_ERROR "COMPILE_MODE must be either debug or optimized" )
    ENDIF()
ENDIF()

# Create the fortran interface file
ADD_DEFINITIONS ( -DCMAKE_CONFIGURED )
IF ( USE_FORTRAN )
    INCLUDE (FortranCInterface)
    FORTRANCINTERFACE_HEADER ( ${AMP_INSTALL_DIR}/include/utils/FC.h MACRO_NAMESPACE "FC_" )
    if( (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") AND (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") )
        set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fno-common")
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fno-common")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-common")
    endif( (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") AND (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin") )
ENDIF()


