INCLUDE(TribitsTplDeclareLibraries)

IF (MSVC)
  # Find the CLAPACK built by CMake on the machine for MSVC if found it will
  # set the BLAS and LAPACK libraries.  NOTE: This the FindCLAPACK module must
  # be called every configure or this does not work!
  FIND_PACKAGE(CLAPACK 3.2.1 NO_MODULE)
  IF (CLAPACK_FOUND AND NOT TPL_BLAS_LIBRARIES)
    ADVANCED_SET(TPL_BLAS_LIBRARIES blas
      CACHE FILEPATH "Set from MSVC CLAPACK specialization")
  ENDIF()
ENDIF()

TRIBITS_TPL_DECLARE_LIBRARIES( BLAS
  REQUIRED_LIBS_NAMES "blas blas_win32"
)

# Add the definitions
SET( USE_EXT_BLAS 1 )
ADD_DEFINITIONS( "-D USE_EXT_BLAS" )
#INCLUDE_DIRECTORIES( ${TPL_BLAS_INCLUDE_DIRS} )
