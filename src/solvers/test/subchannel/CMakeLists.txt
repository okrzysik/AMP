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


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_LIBMESH AND USE_TRILINOS )
    IF( USE_TRILINOS_NOX )
        ADD_TB_TEST( testSubchannelSolveNOX testSubchannelSolveNOX-3 PROCS 4 )
        ADD_TB_TEST( testSubchannelSolveNOX testSubchannelSolveNOX-1 PROCS 8 WEEKLY )
        ADD_TB_TEST( testSubchannelSolutionNOX )
    ENDIF()
ENDIF()
IF( USE_PETSC AND USE_LIBMESH AND USE_TRILINOS )
    ADD_TB_TEST( testSubchannelSolve testSubchannelSolve-3 PROCS 4 )
    ADD_TB_TEST( testSubchannelSolve testSubchannelSolve-1 PROCS 8 WEEKLY )
    ADD_TB_TEST( testSubchannelSolution)
    ADD_TB_TEST( testSubchannelFourEqMMS )
ENDIF()
FINALIZE_TESTBUILDER()



