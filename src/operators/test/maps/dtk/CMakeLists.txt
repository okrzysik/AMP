INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )

COPY_TEST_FILE( input_testDTKEntity )
ADD_AMP_TEST_1_2_4( testDTKEntity )

COPY_TEST_FILE( input_testDTKEntityIterator )
ADD_AMP_TEST_1_2_4( testDTKEntityIterator )

COPY_TEST_FILE( input_testDTKEntitySet )
ADD_AMP_TEST_1_2_4( testDTKEntitySet )

COPY_TEST_FILE( input_testDTKEntityLocalMap )
ADD_AMP_TEST_1_2_4( testDTKEntityLocalMap )

COPY_TEST_FILE( input_testDTKNodalShapeFunction )
ADD_AMP_TEST_1_2_4( testDTKNodalShapeFunction )

COPY_TEST_FILE( input_testDTKAMPMeshManager )
ADD_AMP_TEST_1_2_4( testDTKAMPMeshManager )

COPY_TEST_FILE( input_testDTKAMPField )
ADD_AMP_TEST_1_2_4( testDTKAMPField )

COPY_TEST_FILE( input_testDTKMapOperator )
ADD_AMP_TEST_1_2_4( testDTKMapOperator )

COPY_TEST_FILE( input_testDTKConstruction )
ADD_AMP_PROVISIONAL_TEST( testDTKMapOperatorConstruction )
