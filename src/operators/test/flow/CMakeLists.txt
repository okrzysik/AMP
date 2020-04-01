IF ( AMP_DATA )
    ADD_AMP_TEST_1_2_4 ( testFlowOperator )
    COPY_TEST_FILE( input_testFlowOperator )
    COPY_MESH_FILE ( cube64.with.boundary.labels.e )
ENDIF()

