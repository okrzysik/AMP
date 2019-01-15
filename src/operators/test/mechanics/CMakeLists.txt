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
#     ADD_AMP_TEST_1_2_4( exe_name args )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST ( exe_name args )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST ( exe_name num_cores args )
#
#  To add an input file from the data/ directory
#     COPY_TEST_FILE( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#


# Copy test inputs
COPY_TEST_FILE( input_testLinearMechanics-reset-1 )
COPY_TEST_FILE( input_testNonlinearMechanics-1-normal )
COPY_TEST_FILE( input_testNonlinearMechanics-1-reduced )
COPY_TEST_FILE( input_testNonlinearMechanics_Creep-1-normal )
COPY_TEST_FILE( input_testNonlinearMechanics_Creep-1-reduced )
COPY_TEST_FILE( input_testNonlinearMechanics_Creep_MatLib-1-normal )
COPY_TEST_FILE( input_testNonlinearMechanics_Creep_MatLib-1-reduced )
COPY_TEST_FILE( input_testNonlinearMechanics_NonlinearStrainHardening-1-normal )
COPY_TEST_FILE( input_testNonlinearMechanics_NonlinearStrainHardening-1-reduced )
COPY_TEST_FILE( input_testVonMises_IsotropicKinematicHardening-1-normal )
COPY_TEST_FILE( input_testNonlinearMechanics_Matpro_Fuel_Creep-1-normal )
COPY_TEST_FILE( input_testFrapconCladCreepMaterialModel-1 )
COPY_TEST_FILE( input_testPoroElasticMaterialModel-1 )
COPY_TEST_FILE( input_testLinearMechanics-eigenValues-1 )
COPY_TEST_FILE( input_testNonlinearMechanics-apply-1 )
COPY_TEST_FILE( input_nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1 )
COPY_TEST_FILE( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1 )
COPY_TEST_FILE( input_testNonlinearThermoMechanics-bug-1 )
COPY_TEST_FILE( input_testNonlinearSwellVonMisesMechanics-1 )
COPY_TEST_FILE( input_testNonlinearSwellVonMisesMechanics_MatLib-1 )
COPY_TEST_FILE( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-Oxygen-UO2MSRZC09-1 )
COPY_TEST_FILE( input_testNonlinearMechanics-reset-1 )

COPY_TEST_FILE( input_linearMechanics+ThermalDiffusion-TUI-1 )
COPY_TEST_FILE( input_linearMechanics+ThermalDiffusion-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearMechanics+FickDiffusion-TUI-1 )
COPY_TEST_FILE( input_linearMechanics+FickDiffusion-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearMechanics+ThermalDiffusion+FickDiffusion-TUI-1 )
COPY_TEST_FILE( input_linearMechanics+ThermalDiffusion+FickDiffusion-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Fick+Robin-TUI-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin-TUI-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Fick+Robin-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Thermal+Dirichlet-TUI-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Thermal+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin+Fick+Dirichlet-TUI-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin+Fick+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Thermal+Fick+Dirichlet-TUI-1 )
COPY_TEST_FILE( input_linearBVP-Mechanics+Thermal+Fick+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_FILE(input_testMaterialModels_ElasticDamageThermalStrain-1 )
COPY_TEST_FILE( input_testConstraintsEliminationOperator )


# Copy test data
COPY_MESH_FILE( brick.e )
COPY_MESH_FILE ( cylinder270.e )
COPY_MESH_FILE ( cylinder270.e )


# Add the tests
INITIALIZE_TESTBUILDER()
ADD_TB_TEST_1_2_4( testLinearMechanics_reset )
ADD_TB_TEST_1_2_4( testLinearMechanics_eigenValues )
ADD_TB_TEST_1_2_4( testNonlinearMechanics )
ADD_TB_TEST_1_2_4( testNonlinearMechanics_apply )
ADD_TB_TEST_1_2_4( testNonlinearThermoMechanics )
ADD_TB_TEST_1_2_4( testNonlinearThermoMechanics_bug )
ADD_TB_TEST_1_2_4( testNonlinearThermalOxygenDiffusionWithMechanics )
ADD_TB_TEST_1_2_4( testNonlinearMechanics_reset )
ADD_TB_TEST_1_2_4( testMaterialModels_ElasticDamageThermalStrain )
IF ( USE_EXT_PETSC )
    ADD_TB_TEST_1_2_4 ( testConstraintsEliminationOperator )
ENDIF()
FINALIZE_TESTBUILDER()



