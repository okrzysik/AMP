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

# Copy meshes and input files
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( clad_1x_1pellet.e )
COPY_MESH_FILE ( clad_1x_2pellet.e )
COPY_MESH_FILE ( clad_1x_3pellet.e )
COPY_MESH_FILE ( cube343H27.e )
COPY_MESH_FILE ( pellet_lo_res.e )
COPY_MESH_FILE ( clad_lo_res_4pellet.e )
COPY_MESH_FILE ( multiElementMesh.e )
COPY_TEST_FILE( input_Mesh )
COPY_TEST_FILE( input_SiloIO-1 )
COPY_TEST_FILE( input_SiloIO-2 )
COPY_TEST_FILE( input_SiloIO-3 )
COPY_TEST_FILE( input_Crystal )
COPY_TEST_FILE( input_crystal_mesh )
COPY_TEST_FILE( inputMultiPellet348 )
COPY_TEST_FILE( 17x17_denovo_thermal_1x )
COPY_TEST_FILE( input_AMP_mesh_zoo_1D )
COPY_TEST_FILE( input_AMP_mesh_zoo_2D )
COPY_TEST_FILE( input_AMP_mesh_zoo_3D )

# test the load balancer
ADD_AMP_TEST( test_LoadBalance 1 input_Mesh )
ADD_AMP_TEST( test_LoadBalance 2 input_Mesh )
ADD_AMP_TEST( test_LoadBalance 4 input_Mesh )
ADD_AMP_TEST( test_LoadBalance 4 input_crystal_mesh 4.0 )
ADD_AMP_TEST( test_LoadBalance 512 inputMultiPellet348 )
ADD_AMP_TEST( test_LoadBalance 64000 inputMultiPellet348 )
#ADD_AMP_TEST( test_LoadBalance 64000 17x17_denovo_thermal_1x )

# Add tests
ADD_AMP_TEST( test_BoxMeshHelpers )
ADD_AMP_TEST_1_2_4( test_Mesh )
ADD_AMP_TEST_PARALLEL( test_Mesh 4 input_AMP_mesh_zoo_2D )
ADD_AMP_TEST_PARALLEL( test_Mesh 4 input_AMP_mesh_zoo_3D )
ADD_AMP_PROVISIONAL_TEST( testMeshPerformance )
IF ( USE_TRILINOS_STKMESH )
ENDIF()
IF ( USE_EXT_LIBMESH )
    ADD_AMP_TEST_1_2_4( test_InitializeLibmesh )
ENDIF()
ADD_AMP_TEST_1_2_4( test_StructuredMeshHelpers )

ADD_AMP_TEST( testHex8ElementMapping )
ADD_AMP_TEST( testHex8ElementLinearInterpolation )
ADD_AMP_TEST( testHex8ElementNormalToFaces )
ADD_AMP_TEST( testTriangleProjection )
IF ( USE_EXT_LIBMESH )
    ADD_AMP_TEST_1_2_4 ( test_SiloIO input_SiloIO-1 )
ELSE()
    ADD_AMP_TEST_1_2_4 ( test_SiloIO input_SiloIO-2 )
ENDIF()
IF( USE_EXT_DENDRO AND USE_EXT_PETSC )
    ADD_AMP_TEST_1_2_4( testDendroInterpolation )
    COPY_TEST_FILE( input_testDendroInterpolation  )
    ADD_AMP_TEST_PARALLEL( testDendroSearch 4)
ENDIF()
ADD_AMP_TEST_1_2_4( testMultiMeshSubset )
COPY_TEST_FILE( input_testMultiMeshSubset )
ADD_AMP_TEST_1_2_4( testMultiMeshOwnerRank )
COPY_TEST_FILE( input_testMultiMeshOwnerRank )

