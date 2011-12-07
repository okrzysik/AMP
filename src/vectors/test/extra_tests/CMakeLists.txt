INCLUDE (macros)

IF ( USE_AMP_MESH ) 
    # Copy some mesh files
    COPY_MESH_FILE ( pellet_1x.e )
    COPY_MESH_FILE ( pellet_4x.e )
    # Add a provisional test for testing the performance of makeConsistent
    ADD_AMP_PROVISIONAL_TEST ( testMakeConsistentPerformance )
    # Add a test for testing 64-bit vectors.  WARNING: this test uses a LOT of memory
    COPY_TEST_DATA_FILE ( input-test64bitVectors )
    ADD_AMP_PROVISIONAL_TEST ( test64bitVectors )
ENDIF()


