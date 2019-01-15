############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST ( exe_name args )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL ( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST ( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST ( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#


# Copy test inputs
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-cube )
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-cylinder )
COPY_TEST_FILE( input_testNodeToSegmentConstraintsOperator-pellet )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-1 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-2 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-2bis )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-3 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-4 )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-4bis )
COPY_TEST_FILE( input_testNodeToFaceContactOperator-5 )


# Copy test data
COPY_MESH_FILE( simple_pellet_1x.e )
COPY_MESH_FILE( two_fused_simple_pellets_1x.e )
COPY_MESH_FILE(my_cylinder.e)


# Add the tests

IF( USE_EXT_DENDRO AND USE_EXT_PETSC )
   # R.S.: This needs to be fixed!
    ## does not work in parallel - ADD_AMP_TEST_1_2_4( testNodeToFaceContactOperator )
    ADD_AMP_TEST( testNodeToFaceContactOperator )
    ADD_AMP_TEST_1_2_4( testNodeToSegmentConstraintsOperator )
    ADD_AMP_TEST( testNodeToFaceContactOperator-1 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-2 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-3 )
    ADD_AMP_WEEKLY_TEST( testNodeToFaceContactOperator-4 1 )
    ADD_AMP_WEEKLY_TEST( testNodeToFaceContactOperator-5 1 )
ENDIF()
