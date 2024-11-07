BEGIN_PACKAGE_CONFIG( matrices )

# Add the subdirectories
ADD_PACKAGE_SUBDIRECTORY( testHelpers )

ADD_PACKAGE_SUBDIRECTORY( data )
IF ( USE_HYPRE )
    ADD_PACKAGE_SUBDIRECTORY( data/hypre )
ENDIF()
IF ( USE_TRILINOS )
    ADD_PACKAGE_SUBDIRECTORY( trilinos )
ENDIF()
IF ( USE_PETSC )
    ADD_PACKAGE_SUBDIRECTORY( petsc )
ENDIF()

IF( USE_CUDA )
   ADD_PACKAGE_SUBDIRECTORY( data/cuda )
   SET_DEVICE_COMPILER ( data/cuda )
ENDIF()
IF( USE_HIP )
    ADD_PACKAGE_SUBDIRECTORY( data/hip )
   SET_DEVICE_COMPILER ( data/hip )
ENDIF()

ADD_PACKAGE_SUBDIRECTORY( operations )
ADD_PACKAGE_SUBDIRECTORY( operations/default )
IF( USE_KOKKOS )
    ADD_PACKAGE_SUBDIRECTORY( operations/kokkos )
    SET_DEVICE_COMPILER( operations/kokkos/CSRMatrixOperationsKokkos.cpp )
ENDIF()

SET_DEVICE_COMPILER( CSRMatrix.cpp )

# Set the default matrix type
SET( MATRIX_TYPES "ManagedEpetraMatrix" "NativePetscMatrix" "CSRMatrix" "DenseSerialMatrix" "NULL" )
IF ( NOT DEFINED DEFAULT_MATRIX )
    IF ( USE_TRILINOS )
        SET( DEFAULT_MATRIX "ManagedEpetraMatrix" )
    ELSE()
        SET( DEFAULT_MATRIX "CSRMatrix" )
    ENDIF()
ENDIF()
IF ( NOT "${DEFAULT_MATRIX}" IN_LIST MATRIX_TYPES )
    MESSAGE( FATAL_ERROR "DEFAULT_MATRIX=${DEFAULT_MATRIX} is not a valid matrix type" )
ENDIF()
IF ( "${DEFAULT_MATRIX}" STREQUAL "ManagedEpetraMatrix" AND NOT USE_TRILINOS )
    MESSAGE( FATAL_ERROR "Specified default matrix of ManagedEpetraMatrix, but compiled without Trilinos" )
ENDIF()
IF ( "${DEFAULT_MATRIX}" STREQUAL "NativePetscMatrix" AND NOT USE_PETSC )
    MESSAGE( FATAL_ERROR "Specified default matrix of NativePetscMatrix, but compiled without PETSc" )
ENDIF()
ADD_COMPILE_DEFINITIONS( DEFAULT_MATRIX="${DEFAULT_MATRIX}" )

# Install the matrices target
INSTALL_AMP_TARGET( matrices )

# Add the tests
ADD_SUBDIRECTORY( test )

