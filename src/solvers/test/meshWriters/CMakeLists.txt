############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST( exe_name args )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE( mesh_file )
#


# Copy test inputs



# Copy test data



# Add the tests
ADD_AMP_PROVISIONAL_TEST( writeBox )
ADD_AMP_PROVISIONAL_TEST( write2elementMesh )
ADD_AMP_PROVISIONAL_TEST( write7elementMesh-1 )
ADD_AMP_PROVISIONAL_TEST( write7elementMesh-2 )
ADD_AMP_PROVISIONAL_TEST( write3pointConstrainedBox )
ADD_AMP_PROVISIONAL_TEST( writeDistortedElementMesh )
ADD_AMP_PROVISIONAL_TEST( writeDispValsForPatchTest-2 )
ADD_AMP_PROVISIONAL_TEST( writePlateWithHoleMesh )




