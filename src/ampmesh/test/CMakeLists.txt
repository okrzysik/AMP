INCLUDE (macros)

# Copy meshes and input files
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( clad_1x_1pellet.e )
COPY_MESH_FILE ( clad_1x_2pellet.e )
COPY_MESH_FILE ( clad_1x_3pellet.e )
COPY_MESH_FILE ( cube343H27.e )
COPY_MESH_FILE ( pellet_lo_res.e )
COPY_MESH_FILE ( clad_lo_res_4pellet.e )
COPY_MESH_FILE ( multiElementMesh.e )
COPY_TEST_DATA_FILE ( input_Mesh )
COPY_TEST_DATA_FILE ( input_SiloIO )
COPY_TEST_DATA_FILE ( input_SiloIO-2 )
COPY_TEST_DATA_FILE ( input_Crystal )

# Add tests
ADD_AMP_TEST_1_2_4 ( test_Mesh )
ADD_AMP_TEST_1_2_4 ( test_SiloIO )

IF( USE_DENDRO AND USE_PETSC )
    ADD_AMP_PROVISIONAL_TEST( testDendroSearch )
    COPY_TEST_DATA_FILE ( input_testDendroSearch )
ENDIF()


