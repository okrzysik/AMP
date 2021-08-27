############################################################################
#  Configure tests
############################################################################

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
COPY_TEST_FILE( input_testOxideTimeIntegrator-1 )


# Copy test data
IF ( AMP_DATA )
    COPY_MESH_FILE( unitCubeDirichlet.e )
    COPY_MESH_FILE( simple_linear_thermomech.e )
ENDIF()


# Add the tests
INITIALIZE_TESTBUILDER()
ADD_TB_TEST_1_2_4( testTimeIntegratorParameters )
IF ( USE_EXT_LAPACK )
    ADD_TB_TEST_1_2_4( testOxideTimeIntegrator )
ENDIF()
IF ( USE_EXT_SUNDIALS AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS AND USE_EXT_PETSC AND AMP_DATA )
    ADD_TB_TEST( testIDATimeIntegrator )
    ADD_TB_TEST_1_2_4( testIDA_LinearBVPOperator )
    ADD_TB_TEST_1_2_4( testIDA_NonlinearBVPOperator )
    #ADD_TB_TEST_1_2_4( testIDA_NonlinearColumnOperator ) # Move to provisional until it is running
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearColumnOperator )
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearThermalOxygenDiffusion )
ENDIF()
FINALIZE_TESTBUILDER()

