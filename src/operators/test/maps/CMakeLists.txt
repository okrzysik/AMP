INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )


# Add Jacobian map
ADD_AMP_PROVISIONAL_TEST ( testJacobianMap )
COPY_MESH_FILE ( distortedElementMesh )

# Add node-node tests
ADD_AMP_TEST_1_2_4 ( testNodeToNodeMap inputNodeToNodeMap-1 )           # This input is valid for all processors
#ADD_AMP_TEST_PARALLEL ( testNodeToNodeMap 4 inputNodeToNodeMap-2 )      # This input requires 4 or more processors
COPY_TEST_FILE( inputNodeToNodeMap-1 )
COPY_TEST_FILE( inputNodeToNodeMap-2 )

# Test ScalarZAxis map
ADD_AMP_TEST_1_2_4 ( testScalarZAxisMap )
COPY_MESH_FILE ( TC1-PELL-COARSE.e )
COPY_MESH_FILE ( CLAD-COARSE.e )
COPY_TEST_FILE( inputScalarZAxisMap-1 )
COPY_MESH_FILE(pellet_coarse.e)
COPY_MESH_FILE(pellet_1x.e)
COPY_MESH_FILE(clad_1x_1pellet.e)

# Test ScalarN2GZAxis map
ADD_AMP_TEST_1_2_4 ( testScalarN2GZAxisMap )
COPY_TEST_FILE( inputScalarN2GZAxisMap-1 )

# Test CladToSubchannel map
ADD_AMP_TEST_1_2_4 ( testCladToSubchannelMap )
ADD_AMP_TEST_1_2_4 ( testSubchannelToCladMap )
COPY_TEST_FILE( inputCladToSubchannelMap-1 )
COPY_TEST_FILE( inputSubchannelToCladMap-1 )


# Test map operators
ADD_AMP_TEST_1_2_4 ( testMapOperator )
COPY_TEST_FILE( input_testMapOperator )
COPY_MESH_FILE(halden_1_pellet_9mil.e)
COPY_MESH_FILE(halden_1_clad.e)

# Test strided z axis map
ADD_AMP_PROVISIONAL_TEST ( testStridedZAxisMap )
COPY_TEST_FILE( input_testStridedZAxisMap )

IF( USE_EXT_DTK )
    ADD_AMP_TEST_1_2_4( testMultiDofDTKMapOperator )
    ADD_AMP_TEST_1_2_4( testMultiDofDTKMapOperator-2 )
    COPY_TEST_FILE( input_testMultiDofDTKMapOperator-1 )
    COPY_TEST_FILE( input_testMultiDofDTKMapOperator-2 )
ENDIF()
