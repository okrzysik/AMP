INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )
############################################################################
#  Configure test library
############################################################################
#
#  This sets up a library to be compiled with each test
#  SET ( PACKAGE_TEST_LIB "solver_test" )
#
#  To add source files to the test library
#     ADD_FILES_TO_TEST_LIB ( source_name )

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
#     COPY_TEST_DATA_FILE ( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#

ADD_AMP_TEST ( testMaterial )
ADD_AMP_TEST ( test_UO2_MSRZC_09 )
ADD_AMP_TEST ( testWaterLibrary )
ADD_AMP_TEST ( test_Independent )
ADD_AMP_TEST ( examiner )
ADD_AMP_TEST ( testExtraParams2Eval )
ADD_AMP_TEST ( test_FixedClad )
ADD_AMP_TEST ( test_FixedFuel )

COPY_TEST_DATA_FILE ( inputExaminer )
