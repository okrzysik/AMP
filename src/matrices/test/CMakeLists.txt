INCLUDE (${AMP_SOURCE_DIR}/cmake/macros.cmake)
############################################################################
#  Configure dependencies for the tests
############################################################################
#
#  The TEST_DEP_LIST variable contains the list of AMP packages that these
#  tests depend on

SET ( TEST_DEP_LIST ${TEST_DEP_LIST} matrices vectors )

############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST ( exe_name )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL ( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4 ( exe_name )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST ( exe_name )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST ( exe_name num_cores )
#
#  To add an input file from the data/ directory
#     COPY_TEST_DATA_FILE ( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#

# Copy meshes
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( pellet_lo_res.e )
COPY_MESH_FILE ( clad_1x_1pellet.e )
COPY_MESH_FILE ( clad_1x_2pellet.e )

# Add matrix test
ADD_AMP_TEST_1_2_4( test_Matrix )

