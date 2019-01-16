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
#     ADD_AMP_TEST_1_2_4 ( exe_name args )
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
COPY_TEST_FILE( input_testIDA-LinearBVPOperator-1 )
COPY_TEST_FILE( input_testIDA-NonlinearBVPOperator-1 )
COPY_TEST_FILE( input_testIDA-NonlinearColumnOperator-1 )
COPY_TEST_FILE( input_testIDA-TO-Preprocessor )
COPY_TEST_FILE( input_testIDA-NonlinearOperator-short )
COPY_TEST_FILE( input_testIDA-NonlinearOperator-long )
COPY_TEST_FILE( input_testIDA-NonlinearThermalOxygenDiffusion-1 )
COPY_TEST_FILE( inputIDATimeIntegrator )
COPY_TEST_FILE( input-LinearIDA-1 )
COPY_TEST_FILE( input-NonlinearIDA-1 )
COPY_TEST_FILE( input_testBackwardEulerTimeIntegrator )
COPY_TEST_FILE( input_testOxideTimeIntegrator-1 )


# Copy test data
COPY_MESH_FILE( unitCubeDirichlet.e )
COPY_MESH_FILE( simple_linear_thermomech.e )


# Add the tests
INITIALIZE_TESTBUILDER()
ADD_TB_TEST_1_2_4( testTimeIntegratorParameters )
ADD_TB_TEST_1_2_4( testOxideTimeIntegrator )
IF ( USE_EXT_LIBMESH AND USE_EXT_TRILINOS AND USE_EXT_PETSC )
    ADD_TB_PROVISIONAL_TEST( testBackwardEulerTimeIntegrator )    
ENDIF()
IF ( USE_EXT_SUNDIALS AND USE_EXT_TRILINOS AND USE_EXT_PETSC )
    ADD_TB_TEST( testIDATimeIntegrator )
    ADD_TB_TEST_1_2_4( testIDA_LinearBVPOperator )
    ADD_TB_TEST_1_2_4( testIDA_NonlinearBVPOperator )
    #ADD_TB_TEST_1_2_4( testIDA_NonlinearColumnOperator ) # Move to provisional until it is running
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearColumnOperator )
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearThermalOxygenDiffusion )
ENDIF()
FINALIZE_TESTBUILDER()

