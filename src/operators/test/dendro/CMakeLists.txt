############################################################################
#  Configure tests
############################################################################

# Copy test inputs
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-cube )
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-cylinder )
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-pellet )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-1 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-2 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-2bis )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-3 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-4 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-5 )


# Copy test data
IF ( AMP_DATA )
    COPY_MESH_FILE( simple_pellet_1x.e )
    COPY_MESH_FILE( two_fused_simple_pellets_1x.e )
ENDIF()


# Add the tests
IF( USE_DENDRO AND USE_PETSC )
   # R.S.: This needs to be fixed!
    ## does not work in parallel - ADD_AMP_TEST_1_2_4( testNodeToFaceContactOperator )
    ADD_AMP_TEST( testNodeToFaceContactOperator )
    ADD_AMP_TEST_1_2_4( testNodeToSegmentConstraintsOperator )
    ADD_AMP_TEST( testNodeToFaceContactOperator-1 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-2 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-3 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-4 1 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-5 1 )
ENDIF()
