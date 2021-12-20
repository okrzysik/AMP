# Add a provisional test for testing the performance of makeConsistent
ADD_AMP_TEST_PARALLEL( testMakeConsistentPerformance 4 )
# Add a test for testing 64-bit vectors.  WARNING: this test uses a LOT of memory
IF ( USE_EXT_LIBMESH AND AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_TEST_FILE( input-test64bitVectors )
    ADD_AMP_WEEKLY_TEST( test64bitVectors 1 )
ENDIF()


