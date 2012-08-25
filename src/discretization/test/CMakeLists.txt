INCLUDE (macros)

# Copy meshes and input files
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( clad_1x_1pellet.e )
COPY_MESH_FILE ( cube343H27.e )
COPY_MESH_FILE ( pellet_lo_res.e )
COPY_MESH_FILE ( clad_lo_res_4pellet.e )
#COPY_TEST_DATA_FILE ( input_Mesh )

# Add tests
ADD_AMP_TEST_1_2_4 ( test_Discretization )
IF ( USES_LIBMESH )
    ADD_AMP_TEST_1_2_4 ( test_libmeshElement )
ENDIF()

