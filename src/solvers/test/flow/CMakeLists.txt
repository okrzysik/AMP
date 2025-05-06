############################################################################
#  Configure tests
############################################################################

# Copy test inputs
COPY_TEST_FILE( input_testLinearFlow-1 )
COPY_TEST_FILE( input_testFlowSolution )
COPY_TEST_FILE( input_testThermalRobinFlow-1 )
COPY_TEST_FILE( input_testFlowSolver )
COPY_TEST_FILE( input_testThermalRobinFlow-2 )


# Copy test data
IF ( AMP_DATA )
    COPY_MESH_FILE( cube343H27.e )
    COPY_MESH_FILE( cube343H8.e )
    COPY_MESH_FILE( cube.with.12.boundaries.e )
    COPY_MESH_FILE( cube.with.12.boundaries.2.e )
    COPY_MESH_FILE( halden_1_pellet_9mil.e )
    COPY_MESH_FILE( halden_1_clad.e )
ENDIF()


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_PETSC AND USE_LIBMESH AND USE_TRILINOS AND AMP_DATA )
    ADD_TB_TEST_1_2_4( testFlowSolution )
    ADD_TB_PROVISIONAL_TEST( testLinearFlow )
    ADD_TB_TEST_1_2_4( testThermalRobinFlow )
    ADD_TB_TEST_1_2_4( testFlowSolver )
ENDIF()
FINALIZE_TESTBUILDER()
