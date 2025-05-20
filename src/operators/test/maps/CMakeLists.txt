INITIALIZE_TESTBUILDER()


# Add node-node tests
ADD_TB_TEST_1_2_4( testNodeToNodeMap inputNodeToNodeMap-1 )
COPY_TEST_FILE( inputNodeToNodeMap-1 )


# Test strided z axis map
ADD_TB_TEST( testStridedZAxisMap )
COPY_TEST_FILE( input_testStridedZAxisMap )


# Test maps that use libmesh
IF ( AMP_DATA AND USE_LIBMESH )

    # Add Jacobian map
    ADD_TB_PROVISIONAL_TEST( testJacobianMap )

    # Test ScalarN2GZAxis map
    ADD_TB_TEST_1_2_4( testScalarN2GZAxisMap )
    COPY_TEST_FILE( inputScalarN2GZAxisMap-1 )

    # Test ScalarZAxis map
    ADD_TB_TEST_1_2_4( testScalarZAxisMap )
    COPY_TEST_FILE( inputScalarZAxisMap-1 )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( clad_1x_1pellet.e )

    # Test map operators
    ADD_TB_TEST_1_2_4( testMapOperator )
    COPY_TEST_FILE( input_testMapOperator )

    # Test CladToSubchannel map
    ADD_TB_TEST_1_2_4( testCladToSubchannelMap )
    ADD_TB_TEST_1_2_4( testSubchannelToCladMap )
    COPY_TEST_FILE( inputCladToSubchannelMap-1 )
    COPY_TEST_FILE( inputSubchannelToCladMap-1 )

ENDIF()


# Test DTK map
IF( USE_DTK AND USE_LIBMESH )
    ADD_TB_TEST_1_2_4( testMultiDofDTKMapOperator )
    ADD_TB_TEST_1_2_4( testMultiDofDTKMapOperator_2 )
    COPY_TEST_FILE( input_testMultiDofDTKMapOperator-1 )
    COPY_TEST_FILE( input_testMultiDofDTKMapOperator-2 )
ENDIF()

FINALIZE_TESTBUILDER()
