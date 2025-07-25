# Copy input files
COPY_TEST_FILE( input_testMatVecPerf-1 )
COPY_TEST_FILE( input_testMatVecPerf-2 )
COPY_TEST_FILE( input_testMatVecPerf-3 )
COPY_TEST_FILE( input_testMatCreate-1 )
COPY_TEST_FILE( input_testMatVec-1 )
COPY_TEST_FILE( input_testMatVec-2 )
FILE( COPY "${AMP_SOURCE_DIR}/mesh/test/data/input_AMP_mesh_zoo_3D" DESTINATION . )
FILE( COPY "${AMP_SOURCE_DIR}/mesh/test/data/input_AMP_mesh_zoo_3D_hires" DESTINATION . )


# Copy meshes
IF ( AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
ENDIF()


# Add matrix test
ADD_AMP_TEST_1_2_4( test_Matrix )
ADD_AMP_TEST_1_2_4( test_RectangularMatricesMultiplication )
ADD_AMP_TEST_1_2_4( testMatMultCoarsen )
ADD_AMP_TEST_1_2_4( testMatTranspose GPU )
ADD_AMP_TEST_1_2_4( testMatVecPerf GPU )
ADD_AMP_TEST_1_2_4( testMatMultPerf )
ADD_AMP_TEST_1_2_4( testMatOpDev GPU )
ADD_AMP_TEST_1_2_4( testMatCopyCast GPU )
ADD_AMP_TEST_1_2_4( test_MatMultNonZeros )
ADD_AMP_TEST( test_getRow input_AMP_mesh_zoo_3D )
ADD_AMP_TEST( test_getRow input_AMP_mesh_zoo_3D PROCS 4 )
IF ( USE_DEVICE )
    SET_DEVICE_COMPILER( testMatCopyCast )
    ADD_AMP_TEST_1_2_4( testCompareCSRMatOps GPU )
    SET_DEVICE_COMPILER( testCompareCSRMatOps )
ENDIF()
IF( USE_LIBMESH AND AMP_DATA )
    ADD_AMP_TEST_1_2_4( testMatCreate )
    # This fails for PETSc at present due to a missing interface to getvecs
    #    IF ( USE_TRILINOS OR USE_PETSC ) 
    IF ( USE_TRILINOS ) 
        ADD_AMP_TEST_1_2_4( testMatVec )
    ENDIF()
ENDIF()
IF ( USE_TRILINOS ) 
    ADD_AMP_TEST( test_Epetra_FECrsMatrix_bug PROCS 2 )
ENDIF()

