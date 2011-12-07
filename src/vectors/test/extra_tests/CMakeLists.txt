INCLUDE (macros)

IF ( USE_AMP_MESH ) 
    COPY_MESH_FILE ( pellet_1x.e )
    ADD_AMP_PROVISIONAL_TEST ( testMakeConsistentPerformance )
ENDIF()


