INCLUDE (${AMP_SOURCE_DIR}/cmake/macros.cmake)
############################################################################
#  Configure test library
############################################################################
#
#  This sets up a library to be compiled with each test
#  SET ( PACKAGE_TEST_LIB "solver_test" )
#
#  To add source files to the test library
#     ADD_FILES_TO_TEST_LIB ( source_name )

############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST ( exe_name args )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL ( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4 ( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST ( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST ( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_DATA_FILE ( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#

COPY_MESH_FILE ( cube8.with.boundary.labels.e )

COPY_MESH_FILE ( TC1-PELL-COARSE.e )

COPY_TEST_DATA_FILE ( input_testNonlinear_MatproDensificationMaterialModel-1 )

COPY_TEST_DATA_FILE ( input_nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 (testParameterFactory-1 )
COPY_TEST_DATA_FILE ( input_testParameterFactory-1 )

ADD_AMP_TEST_1_2_4 ( testOperatorParameterUnitTest )


IF ( USE_EXT_LIBMESH )

    ADD_AMP_TEST_1_2_4 ( testElementOperationFactory-1 )
    COPY_TEST_DATA_FILE ( input_testElementOperationFactory-1 )
    COPY_TEST_DATA_FILE ( input_MassLinearElement-create-1 )

    COPY_TEST_DATA_FILE ( input_testSurfaceIteratorBug )
    ADD_AMP_TEST( testSurfaceIteratorBug )

    ADD_AMP_TEST_1_2_4 ( testNeumannVectorCorrection )
    COPY_TEST_DATA_FILE ( input_testNeumannVectorCorrection )

    ADD_AMP_TEST_1_2_4 ( testNonlinearRobin-1 )
    COPY_TEST_DATA_FILE ( input_testNonlinearRobin-1 )
    COPY_MESH_FILE ( cylinder.e )

    ADD_AMP_TEST( testStiffnessTransformation )

    ADD_AMP_TEST( testPolarDecomposition )

    ADD_AMP_PROVISIONAL_TEST( testInverseJacobianBug )

    ADD_AMP_TEST_1_2_4 ( testPowerShapeToVolIntOperator )
    COPY_TEST_DATA_FILE ( input_testPowerShapeToVolIntOperator )

    ADD_AMP_TEST_1_2_4 ( testPowerShape-1 )
    COPY_TEST_DATA_FILE ( input_testPowerShape-1 )
    COPY_TEST_DATA_FILE ( input_testPowerShape-diffusion )

    ADD_AMP_TEST_1_2_4 ( testPowerShape-2 )
    COPY_TEST_DATA_FILE ( input_testPowerShape-2 )

    ADD_AMP_TEST_1_2_4 ( testPowerShape-3 )
    COPY_TEST_DATA_FILE ( input_testPowerShape-3 )

    ADD_AMP_TEST_1_2_4 ( testPowerShape-5 )
    ADD_AMP_TEST_1_2_4 ( testPowerShape-Zr )
    COPY_TEST_DATA_FILE ( input_testPowerShape-5 )

    ADD_AMP_TEST_1_2_4 ( testLinearOperator-1 )
    COPY_TEST_DATA_FILE ( input_testLinearOperator-1 )
    COPY_MESH_FILE ( brick.e )

    ADD_AMP_TEST_1_2_4 ( testLinearRobin-1 )
    COPY_TEST_DATA_FILE ( input_LinearOp-Robin-1 )
    COPY_MESH_FILE ( cube.with.12.boundaries.e )

    ADD_AMP_TEST_1_2_4 ( testVolumeIntegral-1 )
    COPY_TEST_DATA_FILE ( input_testVolumeIntegral-1 )
    COPY_TEST_DATA_FILE ( input_testVolumeIntegral-2 )
    COPY_TEST_DATA_FILE ( input_testVolumeIntegral-3 )
    COPY_MESH_FILE ( cube64.with.boundary.labels.e )

    ADD_AMP_TEST_1_2_4 ( testNonlinearFickSoret-1 )
    COPY_TEST_DATA_FILE ( input_FickSoret-TUI-1 )
    COPY_TEST_DATA_FILE ( input_FickSoret-UO2MSRZC09-1 )

    ADD_AMP_TEST_1_2_4 ( testFickSoretBVP-1 )
    COPY_TEST_DATA_FILE ( input_FickSoret-BVP-TUI-1 )
    COPY_TEST_DATA_FILE ( input_FickSoret-BVP-UO2MSRZC09-1 )

    ADD_AMP_TEST_1_2_4 ( testLinearColumnOperator-1 )
    COPY_TEST_DATA_FILE ( input_testLinearColumnOperator-1 )
    COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Dirichlet-TUI-1 )
    COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Dirichlet-UO2MSRZC09-1 )
    COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Robin-TUI-1 )
    COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Robin-UO2MSRZC09-1 )

    ADD_AMP_TEST_1_2_4 ( testNeutronicsRhs-1 )
    COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-1 )
    COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-2 )
    COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-3 )
    COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-4 )
    COPY_MESH_FILE ( cube64.with.boundary.labels.e )

    ADD_AMP_TEST_1_2_4 ( testNeutronicsRhs-2 )
    COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-db )

    ADD_AMP_TEST_1_2_4 ( testConsistentTangent )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-1-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-1-mesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-1-cookMesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-1-cookMesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-2-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-2-mesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-2-cookMesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-2-cookMesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-3-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-3-mesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-3-cookMesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-3-cookMesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-4-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-4-mesh1-reduced )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-4-cookMesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangent-4-cookMesh1-reduced )
    COPY_MESH_FILE ( mesh1 )
    COPY_MESH_FILE ( cookMesh1 )

    ADD_AMP_TEST_1_2_4 ( testConsistentTangentBVP )
    COPY_TEST_DATA_FILE ( input_testConsistentTangentBVP-1-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangentBVP-2-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangentBVP-3-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangentBVP-4-mesh1-normal )
    COPY_TEST_DATA_FILE ( input_testConsistentTangentBVP-4-mesh1-reduced )

ENDIF()

COPY_MESH_FILE ( distortedElementMesh )


# ADD_AMP_TEST_1_2_4 ( testGapConductanceOperator )

# COPY_TEST_DATA_FILE ( input_testGapConductanceOperator )

COPY_MESH_FILE ( halden_1_pellet_9mil.e )
COPY_MESH_FILE ( halden_1_clad.e )


#NEK-5000 tests
IF( USE_EXT_NEK ) 
  COPY_TEST_DATA_FILE ( SESSION.NAME )
  COPY_TEST_DATA_FILE ( input.h5m )
  COPY_TEST_DATA_FILE ( pipe.rea )
  ADD_AMP_TEST_PARALLEL( testMoabNekPipe 2 ) # This test doesn't work in serial
  ADD_AMP_TEST_PARALLEL( testMoabNekPipe 4 )
  ADD_AMP_TEST_PARALLEL( testMoabNekOperator 2 ) # By extension, neither does this
  ADD_AMP_TEST_PARALLEL( testMoabNekOperator 4 )
ENDIF()

#MOAB tests
IF( USE_EXT_MOAB )
  COPY_TEST_DATA_FILE ( input_testMoabBasedOperator )
  ADD_AMP_TEST_1_2_4( testMoabBasedOperator )
  
  # Can be removed in new mesh design,
  #  I haven't confirmed that I'm ready to delete 
  #  it entirely yet.
  ADD_AMP_PROVISIONAL_TEST( testMoabInterface )
ENDIF()

IF( USE_EXT_DENDRO AND USE_EXT_PETSC )
   # R.S.: This needs to be fixed!
    ## does not work in parallel - ADD_AMP_TEST_1_2_4( testNodeToFaceContactOperator )
    ADD_AMP_TEST( testNodeToFaceContactOperator )
    ADD_AMP_TEST_1_2_4( testNodeToSegmentConstraintsOperator )
    COPY_TEST_DATA_FILE ( input_testNodeToSegmentConstraintsOperator-cube )
    COPY_TEST_DATA_FILE ( input_testNodeToSegmentConstraintsOperator-cylinder )
    COPY_TEST_DATA_FILE ( input_testNodeToSegmentConstraintsOperator-pellet )
    COPY_MESH_FILE( simple_pellet_1x.e )
    COPY_MESH_FILE( two_fused_simple_pellets_1x.e )
    ADD_AMP_TEST( testNodeToFaceContactOperator-1 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-2 )
    ADD_AMP_TEST( testNodeToFaceContactOperator-3 )
    COPY_MESH_FILE(my_cylinder.e)
    ADD_AMP_WEEKLY_TEST( testNodeToFaceContactOperator-4 1 )
    ADD_AMP_WEEKLY_TEST( testNodeToFaceContactOperator-5 1 )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-1 )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-2 )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-2bis )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-3 )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-4 )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-4bis )
    COPY_TEST_DATA_FILE ( input_testNodeToFaceContactOperator-5 )
ENDIF()


IF( USE_EXT_LIBMESH )
    ADD_SUBDIRECTORY( maps )
    ADD_SUBDIRECTORY( diffusion )
    ADD_SUBDIRECTORY( flow )
    ADD_SUBDIRECTORY( subchannel )
    ADD_SUBDIRECTORY( mechanics )
ENDIF()

