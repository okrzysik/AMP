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

ADD_AMP_TEST_1_2_4 ( test_Vector )
ADD_AMP_TEST_1_2_4 ( test_VectorClone )
ADD_AMP_TEST_1_2_4 ( test_VectorCloneView )
ADD_AMP_TEST_1_2_4 ( test_VectorView )
ADD_AMP_TEST_1_2_4 ( test_VectorSelector )
ADD_AMP_TEST_1_2_4 ( testMultiVector )

COPY_TEST_FILE( input_testMultiVector )
COPY_MESH_FILE ( brick.e )

ADD_SUBDIRECTORY ( extra_tests )
