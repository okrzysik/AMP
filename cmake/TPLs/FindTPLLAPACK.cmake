INCLUDE(TribitsTplDeclareLibraries)

IF (MSVC AND NOT TPL_LAPACK_LIBRARIES AND CLAPACK_FOUND)
  ADVANCED_SET(TPL_LAPACK_LIBRARIES lapack
      CACHE FILEPATH "Set from MSVC CLAPACK specialization")
ENDIF()

TRIBITS_TPL_DECLARE_LIBRARIES( LAPACK
  REQUIRED_LIBS_NAMES "lapack lapack_win32")

# Add the definitions
SET( USE_EXT_LAPCK 1 )

# Write the blas/lapack header
SET( BLAS_LAPACK_HEADER ${AMP_INSTALL_DIR}/include/blas_lapack.h )
CONFIGURE_FILE( ${AMP_SOURCE_DIR}/cmake/BlasLapack/fortran_calls.h ${AMP_INSTALL_DIR}/include/fortran_calls.h COPYONLY )
FILE(WRITE ${BLAS_LAPACK_HEADER} "// This is a automatically generated file to include blas/lapack headers\n" )
FILE(APPEND ${BLAS_LAPACK_HEADER} "#ifndef INCLUDE_BLAS_LAPACK\n" )
FILE(APPEND ${BLAS_LAPACK_HEADER} "#define INCLUDE_BLAS_LAPACK\n" )
FILE(APPEND ${BLAS_LAPACK_HEADER} "#include \"${AMP_INSTALL_DIR}/include/fortran_calls.h\"\n" )
FILE(APPEND ${BLAS_LAPACK_HEADER} "#endif\n" )

