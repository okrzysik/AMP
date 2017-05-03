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
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#

# Copy meshes
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( pellet_lo_res.e )
COPY_MESH_FILE ( clad_1x_1pellet.e )
COPY_MESH_FILE ( clad_1x_2pellet.e )
COPY_MESH_FILE ( cube343H27.e )
COPY_MESH_FILE( mpcMesh-1 )

# Add matrix test
ADD_AMP_TEST_1_2_4( test_Matrix )
ADD_AMP_TEST_1_2_4( test_RectangularMatricesMultiplication )

IF( USE_EXT_LIBMESH )
  ADD_AMP_TEST_1_2_4( testMatCreate )
  COPY_TEST_FILE( input_testMatCreate-1 )
  ADD_AMP_TEST_1_2_4( testMatMultiply )
ENDIF()

ADD_SUBDIRECTORY ( extra_tests )

