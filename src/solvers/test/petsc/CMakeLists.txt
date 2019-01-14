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
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermal-cylinder_MATPRO2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermal-cylinder_MATPRO_clad )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermalOxygenDiffusion-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermalOxygenDiffusion-2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1a )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1b )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1c )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1d )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-TUI-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-TUI-2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-3 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-BLog-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermal-cylinder_kIsOne )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( inputPetscKrylovSolver1 )
COPY_TEST_FILE( inputPetscKrylovSolver2 )
COPY_TEST_FILE( inputPetscSNESSolver1 )
COPY_TEST_FILE( inputPetscSNESSolver2 )
COPY_TEST_FILE( inputPetscSNESSolver3 )
COPY_TEST_FILE( inputPetscSNESSolver4 )
COPY_TEST_FILE( petsc_options1 )
COPY_TEST_FILE( petsc_options2 )
COPY_TEST_FILE( petsc_options3 )
COPY_TEST_FILE( petsc_options4 )
COPY_TEST_FILE( petsc_options5 )
COPY_TEST_FILE( petsc_options6 )
COPY_TEST_FILE( petsc_options7 )
COPY_TEST_FILE( petsc_options8 )
COPY_TEST_FILE( petsc_options9 )


# Copy test data
COPY_MESH_FILE( brick.e )
COPY_MESH_FILE( halden_3pel_clad_1x.e )
COPY_MESH_FILE( halden_contact_pellet_9mil_1x.e )
COPY_MESH_FILE( halden_1_pellet_9mil.e )
COPY_MESH_FILE( halden_1_clad.e )
COPY_MESH_FILE( cylinder4K.e )
COPY_MESH_FILE( halden_3pel_clad_1x.e )


# Add the tests that only require petsc
IF ( USE_EXT_TRILINOS )
    ADD_AMP_TEST_1_2_4( testPetscSNESSolver )
    COPY_TEST_FILE( input_testPetscSNESSolver )
ENDIF()


# Add the tests that require libmesh and petsc
IF ( USE_EXT_LIBMESH AND USE_EXT_TRILINOS )
    ADD_AMP_TEST_1_2_4( testPetscSNESSolver-NonlinearThermal-1 )
    ADD_AMP_WEEKLY_TEST( testPetscSNESSolver-NonlinearThermal-2 8 )
    ADD_AMP_WEEKLY_TEST( testPetscSNESSolver-NonlinearFick-1 8 )
    ADD_AMP_WEEKLY_TEST( testPetscSNESSolver-NonlinearFickSoret-1 1)
    ADD_AMP_WEEKLY_TEST( testPetscSNESSolver-NonlinearFickSoret-1 2)
    ADD_AMP_TEST_PARALLEL( testPetscSNESSolver-NonlinearFickSoret-1 4 )
    ADD_AMP_TEST_1_2_4( testPetscSNESSolver-NonlinearFickSoret-2 )
    ADD_AMP_PROVISIONAL_TEST( testPetscSNESSolver-IncompressibleFlow-1 1)
ENDIF()


