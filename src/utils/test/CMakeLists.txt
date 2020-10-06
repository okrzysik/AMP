############################################################################
#  Configure tests
############################################################################

ADD_AMP_TEST( test_compile )
ADD_AMP_TEST_1_2_4( test_AMPManager )
ADD_AMP_TEST_PARALLEL( test_AMPManager 4 2 )
ADD_AMP_TEST_1_2_4( test_AMP_MPI )
ADD_AMP_TEST_1_2_4( test_Utilities )
ADD_AMP_TEST( test_LapackWrappers )
ADD_AMP_TEST_1_2_4( test_InputManager )
ADD_AMP_TEST( test_Database laser_plasma_input.txt )
ADD_AMP_TEST( test_Database2 )
ADD_AMP_TEST( test_Units )

IF( USE_EXT_PETSC AND USE_EXT_TRILINOS )
ADD_AMP_TEST_1_2_4( test_AsciiWriter )
ELSE()
ADD_AMP_TEST( test_AsciiWriter )
ENDIF()

ADD_AMP_TEST_1_2_4( testManufacturedSolution )
ADD_AMP_TEST( test_LinearRegression )
IF( NOT USE_ICC )
ADD_AMP_TEST( testArray )
ENDIF()
ADD_AMP_TEST( test_MathExpr )

COPY_TEST_FILE( input_InputManager )
COPY_TEST_FILE( laser_plasma_input.txt )
COPY_TEST_FILE( input_Database )
COPY_TEST_FILE( input_SAMRAI )

ADD_AMP_PROVISIONAL_TEST( test_crash )

ADD_AMP_TEST( test_kdtree )
ADD_AMP_TEST( test_DelaunayInterpolation )
