IF ( USE_EXT_TRILINOS ) 
    ADD_AMP_TEST_PARALLEL( test_Epetra_FECrsMatrix_bug 2 )
ENDIF()


