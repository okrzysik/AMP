INITIALIZE_TESTBUILDER()

# We need to sort out which tests rely on AMP_DATA
IF ( AMP_DATA )

    # Add Jacobian map
    ADD_TB_PROVISIONAL_TEST( testJacobianMap )
    COPY_MESH_FILE( distortedElementMesh )

    # Add node-node tests
    ADD_TB_TEST_1_2_4 ( testNodeToNodeMap inputNodeToNodeMap-1 )           # This input is valid for all processors
    #ADD_TB_TEST_PARALLEL ( testNodeToNodeMap 4 inputNodeToNodeMap-2 )      # This input requires 4 or more processors
    COPY_TEST_FILE( inputNodeToNodeMap-1 )
    COPY_TEST_FILE( inputNodeToNodeMap-2 )

    # Test ScalarZAxis map
    ADD_TB_TEST_1_2_4( testScalarZAxisMap )
    COPY_MESH_FILE( TC1-PELL-COARSE.e )
    COPY_MESH_FILE( CLAD-COARSE.e )
    COPY_TEST_FILE( inputScalarZAxisMap-1 )
    COPY_MESH_FILE( pellet_coarse.e )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( clad_1x_1pellet.e )

    # Test ScalarN2GZAxis map
    ADD_TB_TEST_1_2_4( testScalarN2GZAxisMap )
    COPY_TEST_FILE( inputScalarN2GZAxisMap-1 )

    # Test CladToSubchannel map
    ADD_TB_TEST_1_2_4( testCladToSubchannelMap )
    ADD_TB_TEST_1_2_4( testSubchannelToCladMap )
    COPY_TEST_FILE( inputCladToSubchannelMap-1 )
    COPY_TEST_FILE( inputSubchannelToCladMap-1 )


    # Test map operators
    ADD_TB_TEST_1_2_4( testMapOperator )
    COPY_TEST_FILE( input_testMapOperator )
    COPY_MESH_FILE(halden_1_pellet_9mil.e)
    COPY_MESH_FILE(halden_1_clad.e)

    # Test strided z axis map
    ADD_TB_PROVISIONAL_TEST ( testStridedZAxisMap )
    COPY_TEST_FILE( input_testStridedZAxisMap )

    # Test DTK map
    IF( USE_EXT_DTK )
        ADD_TB_TEST_1_2_4( testMultiDofDTKMapOperator )
        ADD_TB_TEST_1_2_4( testMultiDofDTKMapOperator_2 )
        COPY_TEST_FILE( input_testMultiDofDTKMapOperator-1 )
        COPY_TEST_FILE( input_testMultiDofDTKMapOperator-2 )
    ENDIF()

ENDIF()

FINALIZE_TESTBUILDER()
