############################################################################
#  Configure tests
############################################################################

# Copy test inputs
COPY_TEST_FILE( input_testAMGuniform)
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-1 )
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-2 )
COPY_TEST_FILE( input_testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_kIsOne )
COPY_TEST_FILE( input_testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_testColumnSolver-LinearThermalAndFick-1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator2 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator3 )
COPY_TEST_FILE( input_testMatrixFreeML-1 )
COPY_TEST_FILE( input_testMeshRefinementDiffusion-1 )
COPY_TEST_FILE( input_testNewFrozenVectorDesign )


# Copy test data
COPY_MESH_FILE( testAMGmesh5 )
COPY_MESH_FILE( simple_pellet_1x.e )
COPY_MESH_FILE( pellet_lo_res.e )
COPY_MESH_FILE( halden_1_pellet_9mil.e )
COPY_MESH_FILE( halden_1_clad.e )
COPY_MESH_FILE( fullMpcMesh-3 )
COPY_MESH_FILE( distortedElementMesh )
COPY_MESH_FILE( cylinder270.e )
COPY_MESH_FILE( cylinder4K.e )
COPY_MESH_FILE( cube4096.e )
COPY_MESH_FILE( cube64.with.boundary.labels.e )
COPY_MESH_FILE( cube.with.12.boundaries.e )
COPY_MESH_FILE( cube.with.12.boundaries.2.e )


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_EXT_TRILINOS AND USE_EXT_LIBMESH )
    ADD_AMP_TEST( testAMGuniform )
ENDIF()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_AMP_TEST( testMatrixFreeML )
    ADD_AMP_PROVISIONAL_TEST( testMeshRefinementDiffusion-1 )
    ADD_AMP_TEST_1_2_4( testLinearRobinBoundaryOperator-1 )
    ADD_AMP_TEST_1_2_4( testColumnSolver-LinearThermalAndFick-1 )
ENDIF()
ADD_AMP_PROVISIONAL_TEST ( testNewFrozenVectorDesign )
FINALIZE_TESTBUILDER()


# Add subdirectories
IF ( USE_EXT_PETSC )
    ADD_SUBDIRECTORY( petsc )
ENDIF()
IF ( USE_EXT_TRILINOS )
    ADD_SUBDIRECTORY( trilinos )
ENDIF()
ADD_SUBDIRECTORY( mechanics )
ADD_SUBDIRECTORY( flow )
ADD_SUBDIRECTORY( subchannel )
ADD_SUBDIRECTORY( meshWriters )

