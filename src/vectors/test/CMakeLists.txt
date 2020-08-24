############################################################################
#  Configure tests
############################################################################

COPY_TEST_FILE( input_testMultiVector )

INITIALIZE_TESTBUILDER()

ADD_TB_TEST_1_2_4( test_Vector )
ADD_TB_TEST_1_2_4( test_VectorClone )
ADD_TB_TEST_1_2_4( test_VectorCloneView )
ADD_TB_TEST_1_2_4( test_VectorView )
ADD_TB_TEST_1_2_4( test_VectorSelector )
ADD_TB_TEST_1_2_4( testMultiVector )

ADD_AMP_TEST( test_VectorPerformance )
ADD_AMP_TEST( test_Scalar )

FINALIZE_TESTBUILDER()


IF ( USE_CUDA )
    ADD_RESOURCE( test_Vector GPU )
    ADD_RESOURCE( test_Vector_2procs GPU )
    ADD_RESOURCE( test_Vector_4procs GPU )
ENDIF()

ADD_SUBDIRECTORY ( extra_tests )
