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
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalOperator-2_HALDEN )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalOperator-2_HALDEN_clad )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalOperator-3 )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalOperator-bar )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalOperator-cylinder )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearFickOperator-bar )
COPY_TEST_FILE( input_testTrilinosMLSolver-LinearThermalRobin )
COPY_TEST_FILE( input_testTrilinosNOX-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( inputTrilinosMLSolver1 )
COPY_TEST_FILE( input_testTrilinosNOXSolver )


# Copy test data
COPY_MESH_FILE( cube27.e )
COPY_MESH_FILE( cube4096.e )
COPY_MESH_FILE( cube64.with.boundary.labels.e )
COPY_MESH_FILE( halden_1_pellet_9mil.e )
COPY_MESH_FILE( halden_1_clad.e )


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_TRILINOS_NOX )
    ADD_TB_TEST_1_2_4( testTrilinosNOXSolver )
ENDIF()
IF( USE_EXT_LIBMESH AND USE_TRILINOS_NOX )
    ADD_TB_TEST_1_2_4( testTrilinosNOX_NonlinearThermal )
ENDIF()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalRobin )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalOperator_bar )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearFickOperator_bar )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalOperator_cylinder )
    ADD_TB_WEEKLY_TEST( testTrilinosMLSolver_LinearThermalOperator_2 8 )
ENDIF()
FINALIZE_TESTBUILDER()

