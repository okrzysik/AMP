############################################################################
#  Configure tests
############################################################################
#
#  Adding tests and input files
#
#  To add a test in serial:
#     ADD_AMP_TEST( exe_name args )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL( exe_name num_cores args )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE( mesh_file )
#


# Copy test inputs
COPY_TEST_FILE( input_testAMGuniform)
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-1 )
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-2 )
COPY_TEST_FILE( input_testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_kIsOne )
COPY_TEST_FILE( input_testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermal-cylinder_MATPRO petsc/data )
COPY_TEST_FILE( input_testPetscSNESSolver-JFNK-ML-NonlinearMechanics-1-normal mechanics/data )
COPY_TEST_FILE( input_testColumnSolver-LinearThermalAndFick-1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator2 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator3 )
COPY_TEST_FILE( input_testMatrixFreeML-1 )
COPY_TEST_FILE( input_testMeshRefinementDiffusion-1 )
COPY_TEST_FILE( input_testNewFrozenVectorDesign )
COPY_TEST_FILE( input_testSolversForUserMatrix-ML )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-TFQMR )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalRobin )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-2 )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-3 )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-bar )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-cylinder )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-shell )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-2_HALDEN )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-2_HALDEN_clad )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearFickOperator-bar )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearElasticityOperator-1 )
COPY_TEST_FILE( input_testSolversForUserMatrix-MueLu )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-TFQMR )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalRobin )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-2 )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-3 )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-bar )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-cylinder )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-shell )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-2_HALDEN )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-2_HALDEN_clad )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearFickOperator-bar )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearElasticityOperator-1 )
COPY_TEST_FILE( input_testSolversForUserMatrix-BoomerAMG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-TFQMR )


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
COPY_MESH_FILE( brick.e )


# Add the tests
INITIALIZE_TESTBUILDER()
IF( USE_EXT_TRILINOS AND USE_EXT_LIBMESH )
    ADD_AMP_TEST( testAMGuniform )
ENDIF()
IF( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_AMP_PROVISIONAL_TEST( testMeshRefinementDiffusion-1 )
    ADD_AMP_TEST( testMatrixFreeML )
    ADD_AMP_TEST_1_2_4( testLinearRobinBoundaryOperator-1 )
    ADD_AMP_TEST_1_2_4( testColumnSolver-LinearThermalAndFick-1 )
    IF ( USE_TRILINOS_MUELU )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearThermalRobin )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearThermalOperator-bar )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearThermalOperator-cylinder )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearThermalOperator-2 )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearFickOperator-bar )
        ADD_AMP_TEST_1_2_4( testTrilinosMueLuSolver-LinearElasticityOperator-1 )
    ENDIF()
    ADD_AMP_TEST( testLinearSolvers-LinearThermalRobin)
    ADD_AMP_TEST_PARALLEL( testLinearSolvers-LinearThermalRobin 2)
ENDIF()
IF ( USE_EXT_LIBMESH AND USE_EXT_HYPRE )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearThermalRobin )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearThermalOperator-bar )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearThermalOperator-cylinder )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearThermalOperator-2 )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearFickOperator-bar )
    ADD_AMP_TEST_1_2_4( testBoomerAMGSolver-LinearElasticityOperator-1 )
    ADD_AMP_TEST_1_2_4( testSolversForUserMatrix )
ENDIF()
IF ( USE_EXT_PETSC AND USE_EXT_LIBMESH AND USE_EXT_HYPRE )
    ADD_AMP_TEST_1_2_4( testPetscSNESSolver-BoomerAMG-NonlinearThermal-1 )
    ADD_AMP_TEST_1_2_4( testPetscSNESSolver-JFNK-BoomerAMG-NonlinearMechanics-1 )
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

