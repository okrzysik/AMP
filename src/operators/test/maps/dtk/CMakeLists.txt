INCLUDE (${AMP_SOURCE_DIR}/cmake/macros.cmake)

COPY_TEST_DATA_FILE( input_testDTKEntity )
ADD_AMP_TEST_1_2_4( testDTKEntity )

COPY_TEST_DATA_FILE( input_testDTKEntityIterator )
ADD_AMP_TEST_1_2_4( testDTKEntityIterator )