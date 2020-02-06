############################################################################
#  Configure tests
############################################################################

# Copy test inputs
COPY_TEST_FILE( input_testSubchannelSolve-1 )
COPY_TEST_FILE( input_testSubchannelSolve-2 )
COPY_TEST_FILE( input_testSubchannelSolve-3 )
COPY_TEST_FILE( input_testSubchannelSolution-1 )
COPY_TEST_FILE( input_testSubchannelSolution-2 )
COPY_TEST_FILE( input_testSubchannelFourEqMMS-1)
COPY_TEST_FILE( input_testSubchannelFourEqMMS-2)
COPY_TEST_FILE( input_testSubchannelSolveNOX-1 )
COPY_TEST_FILE( input_testSubchannelSolveNOX-3 )


# Copy test data


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    IF( USE_TRILINOS_NOX )
        ADD_TB_TEST_PARALLEL( testSubchannelSolveNOX 4 testSubchannelSolveNOX-3 )
        ADD_TB_WEEKLY_TEST( testSubchannelSolveNOX 8 testSubchannelSolveNOX-1 )
        ADD_TB_TEST( testSubchannelSolutionNOX )
    ENDIF()
ENDIF()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_TB_TEST_PARALLEL( testSubchannelSolve 4 testSubchannelSolve-3 )
    ADD_TB_WEEKLY_TEST( testSubchannelSolve 8 testSubchannelSolve-1 )
    ADD_TB_TEST( testSubchannelSolution)
    ADD_TB_TEST( testSubchannelFourEqMMS )
ENDIF()
FINALIZE_TESTBUILDER()



