############################################################################
#  Configure tests
############################################################################

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
IF ( AMP_DATA )
    COPY_MESH_FILE( cube27.e )
    COPY_MESH_FILE( cube4096.e )
    COPY_MESH_FILE( cube64.with.boundary.labels.e )
    COPY_MESH_FILE( halden_1_pellet_9mil.e )
    COPY_MESH_FILE( halden_1_clad.e )
ENDIF()


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_TRILINOS_NOX AND AMP_DATA )
    ADD_TB_TEST_1_2_4( testTrilinosNOXSolver )
ENDIF()
IF( USE_LIBMESH AND USE_TRILINOS_NOX AND AMP_DATA )
    ADD_TB_TEST_1_2_4( testTrilinosNOX_NonlinearThermal )
ENDIF()
IF( USE_PETSC AND USE_LIBMESH AND USE_TRILINOS AND AMP_DATA )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalRobin )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalOperator_bar )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearFickOperator_bar )
    ADD_TB_TEST_1_2_4( testTrilinosMLSolver_LinearThermalOperator_cylinder )
    ADD_TB_WEEKLY_TEST( testTrilinosMLSolver_LinearThermalOperator_2 8 )
ENDIF()
FINALIZE_TESTBUILDER()

