INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )
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

ADD_AMP_TEST_1_2_4 ( testTimeIntegratorParameters )
ADD_AMP_PROVISIONAL_TEST ( testBackwardEulerTimeIntegrator )
ADD_AMP_TEST_1_2_4 ( testOxideTimeIntegrator )

COPY_TEST_FILE( input_testBackwardEulerTimeIntegrator )
COPY_TEST_FILE( input_testOxideTimeIntegrator-1 )

COPY_MESH_FILE ( unitCubeDirichlet.e )
COPY_MESH_FILE ( simple_linear_thermomech.e )
COPY_MESH_FILE ( clad_lo_res_4pellet.e )

# Sundials tests
IF ( USE_EXT_SUNDIALS )
    ADD_AMP_TEST ( testIDATimeIntegrator )
    ADD_AMP_TEST_1_2_4 ( testIDA-LinearBVPOperator-1 )
    ADD_AMP_TEST_1_2_4 ( testIDA-NonlinearBVPOperator-1 )
    #ADD_AMP_TEST_1_2_4 ( testIDA-NonlinearColumnOperator-1 ) # Move to WEEKLY until it is running
    ADD_AMP_PROVISIONAL_TEST ( testIDA-NonlinearColumnOperator-1 1 )
    ADD_AMP_PROVISIONAL_TEST ( testIDA-NonlinearColumnOperator-1 2 )
    ADD_AMP_PROVISIONAL_TEST ( testIDA-NonlinearColumnOperator-1 4 )
    ADD_AMP_PROVISIONAL_TEST ( testIDA-NonlinearThermalOxygenDiffusion-1 )
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
ENDIF()

