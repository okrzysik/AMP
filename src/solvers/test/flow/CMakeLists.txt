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
COPY_TEST_FILE( input_testLinearFlow-1 )
COPY_TEST_FILE( input_testFlowSolution )
COPY_TEST_FILE( input_testThermalRobinFlow-1 )
COPY_TEST_FILE( input_testFlowSolver )
COPY_TEST_FILE( input_testThermalRobinFlow-2 )


# Copy test data
COPY_MESH_FILE( cube343.e )
COPY_MESH_FILE( cube343H27.e )
COPY_MESH_FILE( cube343H8.e )
COPY_MESH_FILE( cube.with.12.boundaries.e )
COPY_MESH_FILE( cube.with.12.boundaries.2.e )
COPY_MESH_FILE( halden_1_pellet_9mil.e )
COPY_MESH_FILE( halden_1_clad.e )


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_TB_TEST_1_2_4( testFlowSolution )
    ADD_TB_PROVISIONAL_TEST( testLinearFlow )
    #ADD_TB_TEST_1_2_4( testThermalRobinFlow_1 )
    ADD_TB_TEST_1_2_4( testThermalRobinFlow_2 )
    ADD_TB_TEST_1_2_4( testFlowSolver )
ENDIF()
FINALIZE_TESTBUILDER()
