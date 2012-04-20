INCLUDE (macros)

############################################################################
#  Configure test library
############################################################################
#
#  This sets up a library to be compiled with each test
#  SET ( PACKAGE_TEST_LIB "operator_test" )
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
#     ADD_AMP_TEST ( exe_name )
#
#  To add a parallel test:
#     ADD_AMP_TEST_PARALLEL ( exe_name , num_cores )
#
#  To add the AMP 1, 2, and 4 processor test
#     ADD_AMP_TEST_1_2_4 ( exe_name )
#
#  To add a build target without placing in the test suite
#     ADD_AMP_PROVISIONAL_TEST ( exe_name )
#
#  To add a parallel test for the weekly regression:
#     ADD_AMP_WEEKLY_TEST ( exe_name num_cores )
#
#  To add an input file from the data/ directory
#     COPY_TEST_DATA_FILE ( input_file )
#
#  To add a mesh from either the data/ directory or the trunk/meshes/ directory
#     COPY_MESH_FILE ( mesh_file )
#

ADD_AMP_PROVISIONAL_TEST(testRemDup)

COPY_TEST_DATA_FILE ( input_testDofMap )
ADD_AMP_PROVISIONAL_TEST(testDofMap)

COPY_TEST_DATA_FILE ( input_testSurfaceIteratorBug )
ADD_AMP_TEST( testSurfaceIteratorBug )

ADD_AMP_PROVISIONAL_TEST( testInverseJacobianBug )

ADD_AMP_TEST( testPolarDecomposition )
ADD_AMP_TEST( testStiffnessTransformation )

COPY_MESH_FILE ( TC1-PELL-COARSE.e )

ADD_AMP_TEST_1_2_4 ( testPowerShapeToVolIntOperator )
COPY_TEST_DATA_FILE ( input_testPowerShapeToVolIntOperator )

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

ADD_AMP_TEST_1_2_4 ( testLinearOperator-1 )
COPY_TEST_DATA_FILE ( input_testLinearOperator-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics-2 )
COPY_TEST_DATA_FILE ( input_linearMechanics-3 )
COPY_TEST_DATA_FILE ( input_linearMechanics-4 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-TUI-Fick-1 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-TUI-Soret-1 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-TUI-Thermal-1 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-UO2MSRZC09-Fick-1 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-UO2MSRZC09-Soret-1 )
COPY_TEST_DATA_FILE ( input_linearDiffusion-UO2MSRZC09-Thermal-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Fick-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Soret-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Thermal-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Fick-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Soret-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Thermal-DirichletMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Fick-RobinMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Soret-RobinMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-TUI-Thermal-RobinMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Fick-RobinMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Soret-RobinMatrixCorrection-1 )
COPY_TEST_DATA_FILE ( input_linearBVPOperator-Diffusion-UO2MSRZC09-Thermal-RobinMatrixCorrection-1 )
COPY_MESH_FILE ( brick.e )

ADD_AMP_TEST_1_2_4 ( testLinearMechanics-reset-1 )
COPY_TEST_DATA_FILE ( input_testLinearMechanics-reset-1 )

ADD_AMP_TEST_1_2_4 ( testLinearMechanics-eigenValues-1 )
COPY_TEST_DATA_FIlE ( input_testLinearMechanics-eigenValues-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearRobin-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearRobin-1 )
COPY_MESH_FILE ( cylinder.e )

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics-1-normal )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics-1-reduced )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_Creep-1-normal )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_Creep-1-reduced )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_Creep_MatLib-1-normal )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_Creep_MatLib-1-reduced )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_NonlinearStrainHardening-1-normal )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_NonlinearStrainHardening-1-reduced )
COPY_TEST_DATA_FILE ( input_testVonMises_IsotropicKinematicHardening-1-normal )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics_Matpro_Fuel_Creep-1-normal )
COPY_TEST_DATA_FILE ( input_testFrapconCladCreepMaterialModel-1 )
COPY_TEST_DATA_FILE ( input_testPoroElasticMaterialModel-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-apply-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics-apply-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermoMechanics-1 )
COPY_TEST_DATA_FILE ( input_nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-UO2MSRZC09-1 )
COPY_MESH_FILE ( cylinder270.e )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermoMechanics-bug-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearThermoMechanics-bug-1 )
COPY_MESH_FILE ( cylinder270.e )

COPY_TEST_DATA_FILE ( input_testNonlinearSwellVonMisesMechanics-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearSwellVonMisesMechanics_MatLib-1 )
COPY_TEST_DATA_FILE ( input_testNonlinear_MatproDensificationMaterialModel-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermalOxygenDiffusion-1 )
COPY_TEST_DATA_FILE ( input_nonlinearBVP-Thermal-Oxygen-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearThermalOxygenDiffusionWithMechanics-1 )
COPY_TEST_DATA_FILE ( input_nonlinearBVP-Mechanics-ThermalStrain-Thermal-Oxygen-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearMechanics-reset-1 )
COPY_TEST_DATA_FILE ( input_testNonlinearMechanics-reset-1 )

ADD_AMP_TEST_1_2_4 ( testElementOperationFactory-1 )
COPY_TEST_DATA_FILE ( input_testElementOperationFactory-1 )
COPY_TEST_DATA_FILE ( input_MechanicsLinearElement-create-1 )
COPY_TEST_DATA_FILE ( input_MechanicsNonlinearElement-create-1 )
COPY_TEST_DATA_FILE ( input_DiffusionLinearElement-create-1 )
COPY_TEST_DATA_FILE ( input_DiffusionNonlinearElement-create-1 )
COPY_TEST_DATA_FILE ( input_MassLinearElement-create-1 )

ADD_AMP_TEST_1_2_4 (testParameterFactory-1 )
COPY_TEST_DATA_FILE ( input_testParameterFactory-1 )

ADD_AMP_TEST_1_2_4 ( testOperatorParameterUnitTest )

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

ADD_AMP_TEST_1_2_4 ( testLinearDiffusion-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-1  )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Fick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Soret-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Thermal-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Fick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Soret-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-TensorFick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-CylindricalFick-1 )

ADD_AMP_TEST_1_2_4 ( testLinearDiffusion-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-2  )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Fick-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Soret-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Thermal-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Fick-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Soret-2 )


ADD_AMP_TEST_1_2_4 ( testLinearRobin-1 )
COPY_TEST_DATA_FILE ( input_LinearOp-Robin-1 )
COPY_MESH_FILE ( cube.with.12.boundaries.e )

ADD_AMP_TEST_1_2_4 ( testFlowOperator )
COPY_TEST_DATA_FILE ( input_testFlowOperator )
COPY_MESH_FILE ( cube64.with.boundary.labels.e )

ADD_AMP_TEST_1_2_4 ( testVolumeIntegral-1 )
COPY_TEST_DATA_FILE ( input_testVolumeIntegral-1 )
COPY_TEST_DATA_FILE ( input_testVolumeIntegral-2 )
COPY_TEST_DATA_FILE ( input_testVolumeIntegral-3 )
COPY_MESH_FILE ( cube64.with.boundary.labels.e )

ADD_AMP_TEST_1_2_4 ( testNonlinearDiffusion-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-1  )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Fick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Soret-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Thermal-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Fick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Soret-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Fick-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Soret-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Thermal-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Fick-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Soret-ActiveTemperatureAndConcentration-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-TensorFick-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-CylindricalFick-1 )

ADD_AMP_TEST_1_2_4 ( testNonlinearDiffusion-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-1  )

ADD_AMP_TEST_1_2_4 ( testNonlinearFickSoret-1 )
COPY_TEST_DATA_FILE ( input_FickSoret-TUI-1 )
COPY_TEST_DATA_FILE ( input_FickSoret-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testDiffusionBVP-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Thermal-2  )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Fick-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-TUI-Soret-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Thermal-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Fick-2 )
COPY_TEST_DATA_FILE ( input_Diffusion-UO2MSRZC09-Soret-2 )
COPY_MESH_FILE ( cube64.with.boundary.labels.e )


ADD_AMP_TEST_1_2_4 ( testFickSoretBVP-1 )
COPY_TEST_DATA_FILE ( input_FickSoret-BVP-TUI-1 )
COPY_TEST_DATA_FILE ( input_FickSoret-BVP-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testDiffusionManufacturedSolution-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-Fick-OxMSRZC09-MMS-1 )
COPY_TEST_DATA_FILE ( input_Diffusion-Fick-TUI-MMS-1 )
COPY_MESH_FILE ( cube4096.e )

ADD_AMP_TEST_1_2_4 ( testTensorDiffusionManufacturedSolution-1 )
COPY_TEST_DATA_FILE ( input_TensorDiffusion-Fick-MMS-1 )
COPY_TEST_DATA_FILE ( input_TensorDiffusion-Fick-MMS-2 )
COPY_MESH_FILE ( cube4096.e )
COPY_MESH_FILE ( shell3150mesh.e )

ADD_AMP_TEST_1_2_4 ( testLinearColumnOperator-1 )
COPY_TEST_DATA_FILE ( input_testLinearColumnOperator-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+ThermalDiffusion-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+ThermalDiffusion-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+FickDiffusion-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+FickDiffusion-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearThermalDiffusion+FickDiffusion-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearThermalDiffusion+FickDiffusion-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+ThermalDiffusion+FickDiffusion-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearMechanics+ThermalDiffusion+FickDiffusion-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Dirichlet-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Robin-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Thermal+Fick+Robin-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Fick+Robin-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Fick+Robin-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Thermal+Dirichlet-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Thermal+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin+Fick+Dirichlet-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Dirichlet+Thermal+Robin+Fick+Dirichlet-UO2MSRZC09-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Thermal+Fick+Dirichlet-TUI-1 )
COPY_TEST_DATA_FILE ( input_linearBVP-Mechanics+Thermal+Fick+Dirichlet-UO2MSRZC09-1 )

ADD_AMP_TEST_1_2_4 ( testNeutronicsRhs-1 )
COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-1 )
COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-2 )
COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-3 )
COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-4 )
COPY_MESH_FILE ( cube64.with.boundary.labels.e )

ADD_AMP_TEST_1_2_4 ( testNeutronicsRhs-2 )
COPY_TEST_DATA_FILE ( input_testNeutronicsRhs-db )

ADD_AMP_TEST_1_2_4( testMultiVector )
COPY_TEST_DATA_FILE ( input_testMultiVector )

ADD_AMP_PROVISIONAL_TEST ( testJacobianMap )
COPY_MESH_FILE ( distortedElementMesh )

ADD_AMP_PROVISIONAL_TEST ( testDiffusionMMSForward )
ADD_AMP_TEST_1_2_4 ( testMaterialModels_ElasticDamageThermalStrain-1 )
COPY_TEST_DATA_FILE (input_testMaterialModels_ElasticDamageThermalStrain-1 )

# ADD_AMP_TEST_1_2_4 ( testGapConductanceOperator )

# COPY_TEST_DATA_FILE ( input_testGapConductanceOperator )

COPY_MESH_FILE ( halden_1_pellet_9mil.e )
COPY_MESH_FILE ( halden_1_clad.e )

# Add node-node tests
ADD_AMP_TEST_1_2_4 ( testNodeToNodeMap inputNodeToNodeMap-1 )           # This input is valid for all processors
ADD_AMP_TEST_PARALLEL ( testNodeToNodeMap 4 inputNodeToNodeMap-2 )      # This input requires 4 or more processors
COPY_TEST_DATA_FILE ( inputNodeToNodeMap-1 )
COPY_TEST_DATA_FILE ( inputNodeToNodeMap-2 )

# Test ScalarZAxis map
ADD_AMP_TEST_1_2_4 ( testScalarZAxisMap )
COPY_MESH_FILE ( TC1-PELL-COARSE.e )
COPY_MESH_FILE ( CLAD-COARSE.e )
COPY_TEST_DATA_FILE ( inputScalarZAxisMap-1 )
COPY_MESH_FILE(pellet_coarse.e)
COPY_MESH_FILE(pellet_1x.e)
COPY_MESH_FILE(clad_1x_1pellet.e)

#NEK-5000 tests
IF( USE_NEK ) 
  COPY_TEST_DATA_FILE ( SESSION.NAME )
  COPY_TEST_DATA_FILE ( input.h5m )
  COPY_TEST_DATA_FILE ( pipe.rea )
  ADD_AMP_TEST_PARALLEL( testNekPipe 2 ) # This test doesn't work in serial
  ADD_AMP_TEST_PARALLEL( testNekPipe 4 )
  ADD_AMP_TEST_PARALLEL( testNekOperator 2 )
  ADD_AMP_TEST_PARALLEL( testNekOperator 4 )
ENDIF()

#MOAB tests
IF( USE_MOAB )
  ADD_AMP_TEST_1_2_4( testMoabInterface )
ENDIF()
