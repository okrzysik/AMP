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
        ADD_TB_TEST( testSubchannelSolutionNOX )
        ADD_TB_TEST( testSubchannelSolveNOX testSubchannelSolveNOX-3 PROCS 4 )
        IF ( NOT ( CMAKE_BUILD_TYPE STREQUAL "Debug" ) )
            ADD_TB_TEST( testSubchannelSolveNOX testSubchannelSolveNOX-1 PROCS 8 )
        ENDIF()
    ENDIF()
    IF( USE_PETSC )
        ADD_TB_TEST( testSubchannelSolution)
        ADD_TB_TEST( testSubchannelFourEqMMS )
        ADD_TB_TEST( testSubchannelSolve testSubchannelSolve-3 PROCS 4 )
        IF ( NOT ( CMAKE_BUILD_TYPE STREQUAL "Debug" ) )
            ADD_TB_TEST( testSubchannelSolve testSubchannelSolve-1 PROCS 8 )
        ENDIF()
    ENDIF()
ENDIF()
FINALIZE_TESTBUILDER()



