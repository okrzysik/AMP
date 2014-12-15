INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE (${AMP_SOURCE_DIR}/cmake/libraries.cmake)
SET_AMP_PACKAGE_FLAGS()

IF ( USE_EXT_TRILINOS ) 
    ADD_AMP_TEST_PARALLEL( test_Epetra_FECrsMatrix_bug 2 )
ENDIF()


