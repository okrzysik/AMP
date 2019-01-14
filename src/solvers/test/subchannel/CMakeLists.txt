############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST( exe_name args )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE( mesh_file )
#


# Copy test inputs
COPY_TEST_FILE( input_testSubchannelSolve-1 )
COPY_TEST_FILE( input_testSubchannelSolve-2 )
COPY_TEST_FILE( input_testSubchannelSolve-3 )
COPY_TEST_FILE( input_testSubchannelSolution-1 )
COPY_TEST_FILE( input_testSubchannelSolution-2 )
COPY_TEST_FILE( input_testSubchannelFourEqMMS-1)
COPY_TEST_FILE( input_testSubchannelFourEqMMS-2)


# Copy test data


# Add the tests
IF( USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    IF( USE_TRILINOS_NOX )
        ADD_AMP_PROVISIONAL_TEST(testSubchannelSolveNOX 8 testSubchannelSolveNOX-1 )
        COPY_TEST_FILE( input_testSubchannelSolveNOX-1 )
        ADD_AMP_TEST(testSubchannelSolutionNOX)
    ENDIF()
ENDIF()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_AMP_TEST_PARALLEL(testSubchannelSolve 4 testSubchannelSolve-3 )
    ADD_AMP_WEEKLY_TEST(testSubchannelSolve 8 testSubchannelSolve-1 )
    ADD_AMP_TEST(testSubchannelSolution)
    ADD_AMP_PROVISIONAL_TEST( testSubchannelFourEqMMS )
ENDIF()
