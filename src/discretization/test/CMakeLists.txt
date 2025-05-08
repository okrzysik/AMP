############################################################################
#  Configure tests
############################################################################

ADD_AMP_TEST_1_2_4( test_Discretization )
ADD_AMP_TEST_1_2_4( test_Neighbors )
IF ( USE_LIBMESH )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( pellet_lo_res.e )
    COPY_MESH_FILE( clad_1x_1pellet.e )
    ADD_AMP_TEST_1_2_4( test_libmeshElement )
ENDIF()
