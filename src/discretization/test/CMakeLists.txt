############################################################################
#  Configure tests
############################################################################

# Copy meshes
IF ( AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( clad_1x_1pellet.e )
    COPY_MESH_FILE( cube343H27.e )
    COPY_MESH_FILE( pellet_lo_res.e )
    COPY_MESH_FILE( clad_lo_res_4pellet.e )
ENDIF()

# Add tests
ADD_AMP_TEST_1_2_4( test_Discretization )
IF ( USE_EXT_LIBMESH )
    ADD_AMP_TEST_1_2_4( test_libmeshElement )
ENDIF()
