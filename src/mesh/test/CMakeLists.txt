############################################################################
#  Configure tests
############################################################################


# Set default number of processors for tests
SET( MESH_PROCS 1 )
IF ( USE_MPI )
    SET( MESH_PROCS 1 )
ENDIF()


# Copy input files
COPY_TEST_FILE( input_libMesh )
COPY_TEST_FILE( input_crystal )
COPY_TEST_FILE( input_crystal_mesh )
COPY_TEST_FILE( input_crystal_triangle )
COPY_TEST_FILE( inputMultiPellet348 )
COPY_TEST_FILE( 17x17_denovo_thermal_1x )
COPY_TEST_FILE( input_AMP_mesh_zoo_1D )
COPY_TEST_FILE( input_AMP_mesh_zoo_2D )
COPY_TEST_FILE( input_AMP_mesh_zoo_3D )
COPY_TEST_FILE( input_AMP_mesh_zoo_3D_hires )
COPY_TEST_FILE( input_circle_frustrum )
COPY_TEST_FILE( input_compare )
COPY_TEST_FILE( input_compare_libmesh )
COPY_TEST_FILE( input_stl )
COPY_TEST_FILE( input_crystal_circular_tri )
COPY_TEST_FILE( CRYSTAL_FINAL_2mm.STL )
COPY_TEST_FILE( input_testMultiMeshOwnerRank )
COPY_TEST_FILE( input_testMultiMeshSubset )
COPY_TEST_FILE( input_LogicalMesh )


# Copy meshes
IF ( AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( Utah_teapot.stl )
ENDIF()


# test the load balancer
ADD_AMP_TEST( test_LoadBalance 8 inputMultiPellet348 NO_RESOURCES )
ADD_AMP_TEST( test_LoadBalance 4 input_crystal_mesh 4.0 NO_RESOURCES )
ADD_AMP_TEST( test_LoadBalance 512 inputMultiPellet348 NO_RESOURCES )
ADD_AMP_TEST( test_LoadBalance 64000 inputMultiPellet348 NO_RESOURCES )
ADD_AMP_TEST( test_LoadBalance 1000000 17x17_denovo_thermal_1x NO_RESOURCES )
IF ( USE_LIBMESH AND AMP_DATA )
    ADD_AMP_TEST( test_LoadBalance 4 input_libMesh )
ENDIF()


# Add mesh tests
ADD_AMP_TEST_1_2_4( test_Mesh )
ADD_AMP_TEST( test_Mesh input_AMP_mesh_zoo_1D )
ADD_AMP_TEST( test_Mesh input_AMP_mesh_zoo_2D PROCS ${MESH_PROCS} )
ADD_AMP_TEST( test_Mesh input_AMP_mesh_zoo_3D PROCS ${MESH_PROCS} )
ADD_AMP_TEST( test_Mesh input_LogicalMesh )
ADD_AMP_TEST( test_Mesh input_crystal_circular_tri )  # triangle mesh (1 rank only)
IF ( AMP_DATA )
    ADD_AMP_TEST( test_Mesh input_stl )  # triangle mesh (1 rank only)
    # ADD_AMP_TEST( test_Mesh input_crystal_triangle )  # triangle mesh (1 rank only)
ENDIF()
ADD_AMP_TEST( testMeshPerformance )
ADD_AMP_TEST( compare_Mesh input_compare )
ADD_AMP_TEST( view_STL CRYSTAL_FINAL_2mm.STL )
ADD_AMP_TEST_1_2_4( test_StructuredMeshHelpers )
ADD_AMP_TEST_1_2_4( testMultiMeshSubset )
ADD_AMP_TEST_1_2_4( testMultiMeshOwnerRank )
ADD_AMP_TEST( test_BoxMeshIndex )


# Write some old meshes
ADD_AMP_TEST( writeMesh )

# Misc tests
ADD_AMP_TEST( testTriangleProjection )
ADD_AMP_TEST( testHex8Element )

