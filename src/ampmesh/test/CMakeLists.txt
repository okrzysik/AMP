INCLUDE (macros)

# Copy meshes
COPY_MESH_FILE ( pellet_1x.e )
COPY_MESH_FILE ( cube343H27.e )

# Add tests
ADD_AMP_TEST_1_2_4 ( test_Mesh )


