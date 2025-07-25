############################################################################
#  Configure tests
############################################################################

# Copy test inputs
COPY_TEST_FILE( input_testIDA-LinearBVPOperator-1 )
COPY_TEST_FILE( input_testIDA-NonlinearBVPOperator-1 )
COPY_TEST_FILE( input_testIDA-NonlinearColumnOperator-1 )
COPY_TEST_FILE( input_testIDA-TO-Preprocessor )
COPY_TEST_FILE( input_testIDA-NonlinearThermalOxygenDiffusion-1 )
COPY_TEST_FILE( inputIDATimeIntegrator )
COPY_TEST_FILE( input_testOxideTimeIntegrator-1 )
COPY_TEST_FILE( input_testImplicitIntegrationWithUserOperator )

# Add the tests
ADD_AMP_TEST( testImplicitIntegrationWithUserOperator )

INITIALIZE_TESTBUILDER()
ADD_TB_TEST( testSimpleTimeIntegration )
IF ( USE_HDF5 )
    ADD_TB_TEST_1_2_4( testSimpleRestartTimeIntegration )
ENDIF()
ADD_TB_TEST_1_2_4( testTimeIntegratorParameters )
IF ( USE_LAPACK_WRAPPERS )
    IF ( USE_LAPACK )
        ADD_TB_TEST_1_2_4( testOxideTimeIntegrator )
    ELSE()
        ADD_TB_PROVISIONAL_TEST( testOxideTimeIntegrator )
    ENDIF()
ENDIF()
IF ( USE_SUNDIALS AND USE_LIBMESH AND USE_TRILINOS AND USE_PETSC AND AMP_DATA )
    ADD_TB_TEST( testIDATimeIntegrator )
    ADD_TB_TEST_1_2_4( testIDA_LinearBVPOperator )
    ADD_TB_TEST_1_2_4( testIDA_NonlinearBVPOperator )
    #ADD_TB_TEST_1_2_4( testIDA_NonlinearColumnOperator ) # Move to provisional until it is running
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearColumnOperator )
    ADD_TB_PROVISIONAL_TEST( testIDA_NonlinearThermalOxygenDiffusion )
ENDIF()
FINALIZE_TESTBUILDER()

