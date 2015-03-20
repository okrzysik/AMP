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

ADD_AMP_TEST( test_compile )
ADD_AMP_TEST_1_2_4( test_AMPManager )
ADD_AMP_TEST_PARALLEL( test_AMPManager 4 2 )
ADD_AMP_TEST_1_2_4( test_AMP_MPI )
ADD_AMP_TEST_1_2_4( test_Utilities )
ADD_AMP_TEST( test_LapackWrappers )
ADD_AMP_TEST_1_2_4( test_InputManager-1 )
ADD_AMP_TEST_1_2_4( test_Database-1 )
ADD_AMP_TEST_1_2_4( test_AsciiWriter )
ADD_AMP_TEST_1_2_4( testManufacturedSolution )
ADD_AMP_TEST( test_LinearRegression )

COPY_TEST_FILE( input_InputManager-1 )
COPY_TEST_FILE( input_Database-1 )




