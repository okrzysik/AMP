INCLUDE (macros)

ADD_AMP_TEST_1_2_4 ( test_Vector )
ADD_AMP_TEST_1_2_4 ( test_VectorClone )
ADD_AMP_TEST_1_2_4 ( test_VectorCloneView )
ADD_AMP_TEST_1_2_4 ( test_VectorView )
ADD_AMP_TEST_1_2_4 ( test_VectorSelector )
ADD_AMP_TEST_1_2_4( testMultiVector )

COPY_TEST_DATA_FILE ( input_testMultiVector )
COPY_MESH_FILE ( brick.e )

ADD_SUBDIRECTORY ( extra_tests )
