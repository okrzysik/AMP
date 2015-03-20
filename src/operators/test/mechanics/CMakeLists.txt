INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )

COPY_MESH_FILE( brick.e )

ADD_AMP_TEST_1_2_4 ( testLinearMechanics-reset-1 )
COPY_TEST_FILE( input_testLinearMechanics-reset-1 )

ADD_AMP_TEST_1_2_4 ( testLinearMechanics-eigenValues-1 )
COPY_TEST_FILE( input_testLinearMechanics-eigenValues-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-1 )
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

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-apply-1 )
COPY_TEST_FILE( input_testNonlinearMechanics-apply-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermoMechanics-1 )
COPY_TEST_FILE( input_nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1 )
COPY_TEST_FILE( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1 )
COPY_MESH_FILE ( cylinder270.e )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermoMechanics-bug-1 )
COPY_TEST_FILE( input_testNonlinearThermoMechanics-bug-1 )
COPY_MESH_FILE ( cylinder270.e )

COPY_TEST_FILE( input_testNonlinearSwellVonMisesMechanics-1 )
COPY_TEST_FILE( input_testNonlinearSwellVonMisesMechanics_MatLib-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermalOxygenDiffusionWithMechanics-1 )
COPY_TEST_FILE( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-Oxygen-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-reset-1 )
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

ADD_AMP_TEST_1_2_4 ( testMaterialModels_ElasticDamageThermalStrain-1 )
COPY_TEST_FILE(input_testMaterialModels_ElasticDamageThermalStrain-1 )

IF ( USE_EXT_PETSC )
    ADD_AMP_TEST_1_2_4 ( testConstraintsEliminationOperator )
    COPY_TEST_FILE( input_testConstraintsEliminationOperator )
ENDIF()




