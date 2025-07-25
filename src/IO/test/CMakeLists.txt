# Copy test inputs
COPY_TEST_FILE( input_SiloIO-1        data )
COPY_TEST_FILE( input_SiloIO-2        data )
COPY_TEST_FILE( input_crystal         ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_stl             ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_AMP_mesh_zoo_1D ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_AMP_mesh_zoo_2D ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_AMP_mesh_zoo_3D ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_crystal_circular_tri ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( CRYSTAL_FINAL_2mm.STL ${AMP_SOURCE_DIR}/mesh/test/data )
COPY_TEST_FILE( input_RK12-1 ${AMP_SOURCE_DIR}/time_integrators/test/data )

# Copy meshes
IF ( AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( Utah_teapot.stl )
ENDIF()


# Add Writer tests
ADD_AMP_TEST( test_Writer input_crystal )
ADD_AMP_TEST( test_Writer input_AMP_mesh_zoo_2D )
ADD_AMP_TEST( test_Writer input_AMP_mesh_zoo_3D )
ADD_AMP_TEST( test_Writer input_crystal PROCS 2 )
IF ( AMP_DATA )
    ADD_AMP_TEST( test_Writer input_stl )
ENDIF()
IF ( USE_LIBMESH AND AMP_DATA )
    ADD_AMP_TEST_1_2_4( test_Writer input_SiloIO-1 )
ENDIF()
ADD_AMP_TEST_1_2_4( test_Writer input_SiloIO-2 )


# Add HDF5 tests
IF ( USE_HDF5 )
    ADD_AMP_EXECUTABLE( view_HDF5 )
    ADD_AMP_TEST( test_HDF5 )
    ADD_AMP_TEST( view_HDF5 test_HDF5.hdf5 DEPENDS test_HDF5 )
    ADD_AMP_TEST( test_XDMF )
    ADD_AMP_TEST( test_XDMF PROCS 2 )
    ADD_AMP_TEST( test_RestartManager input_AMP_mesh_zoo_1D )
    ADD_AMP_TEST( test_RestartManager input_AMP_mesh_zoo_2D )
    ADD_AMP_TEST( test_RestartManager input_AMP_mesh_zoo_3D )
    ADD_AMP_TEST( test_RestartManagerComms )
    ADD_AMP_TEST( test_RestartManagerComms PROCS 4 )
    ADD_AMP_TEST( test_RestartManager input_RK12-1 )
ENDIF()



