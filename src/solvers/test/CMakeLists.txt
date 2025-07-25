###########################################################################
############################################################################
#  Configure tests

# Copy solver parameter files
COPY_TEST_FILE( CGParameters )
COPY_TEST_FILE( GMRESParameters )
COPY_TEST_FILE( FGMRESParameters )
COPY_TEST_FILE( BiCGSTABParameters )
COPY_TEST_FILE( TFQMRParameters )
COPY_TEST_FILE( GMRESRParameters )
COPY_TEST_FILE( PetscFGMRESParameters )
COPY_TEST_FILE( PetscCGParameters )
COPY_TEST_FILE( PetscBiCGSTABParameters )
COPY_TEST_FILE( PC-BoomerAMGParameters )
COPY_TEST_FILE( PC-MLParameters )
COPY_TEST_FILE( PC-MueLuParameters )
COPY_TEST_FILE( PC-DiagParameters )

# Copy test inputs
COPY_TEST_FILE( input_LinearElasticityOperator-1 )

COPY_TEST_FILE( input_LinearFickOperator-bar )

COPY_TEST_FILE( input_LinearThermalOperator-bar )
COPY_TEST_FILE( input_LinearThermalRobinOperator )
COPY_TEST_FILE( input_LinearThermalRobinOperator-CylMesh )

COPY_TEST_FILE( input_NonlinearFick-cylinder-1a )
COPY_TEST_FILE( input_NonlinearFick-cylinder-1b )
COPY_TEST_FILE( input_NonlinearFick-cylinder-1c )
COPY_TEST_FILE( input_NonlinearFick-cylinder-1d )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-TUI-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-TUI-2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-1 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-2 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-TUI-3 )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFickSoret-cylinder-OxMSRZC09-1 )

COPY_TEST_FILE( input_PetscSNESSolver)
COPY_TEST_FILE( input_TrilinosNOXSolver)

COPY_TEST_FILE( input_testAMGuniform)
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-1 )
COPY_TEST_FILE( input_testLinearRobinBoundaryOperator-2 )
COPY_TEST_FILE( input_testNonlinearKrylovAccelerator-NonlinearThermal-cylinder_kIsOne )
COPY_TEST_FILE( input_NKA-BoomerAMG-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_NKA-TrilinosML-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_PetscSNESSolver-NonlinearThermal-cylinder_MATPRO petsc/data )
COPY_TEST_FILE( input_PetscSNESSolver-NonlinearThermal-cylinder_MATPRO2 petsc/data )
COPY_TEST_FILE( input_PetscSNESSolver-BoomerAMG-NonlinearThermal-cylinder_MATPRO2 )
COPY_TEST_FILE( input_PetscSNESSolver-BoomerAMG-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_PetscSNESSolver-TrilinosML-NonlinearThermal-cylinder_MATPRO )
COPY_TEST_FILE( input_TrilinosNOX-TrilinosML-NonlinearThermal-cylinder_MATPRO )

COPY_TEST_FILE( input_testPetscSNESSolver-JFNK-ML-NonlinearMechanics-1-normal mechanics/data )
COPY_TEST_FILE( input_testPetscSNESSolver-JFNK-BoomerAMG-NonlinearMechanics-1-normal mechanics/data )

COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1a )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1b )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1c )
COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearFick-cylinder-1d )

COPY_TEST_FILE( input_NKASolver-ML-NonlinearMechanics-1-normal mechanics/data )
COPY_TEST_FILE( input_NKASolver-BoomerAMG-NonlinearMechanics-1-normal mechanics/data )

COPY_TEST_FILE( input_testColumnSolver-LinearThermalAndFick-1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator1 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator2 )
COPY_TEST_FILE( inputNonlinearKrylovAccelerator3 )
COPY_TEST_FILE( input_testMatrixFreeML-1 )
COPY_TEST_FILE( input_testMeshRefinementDiffusion-1 )
COPY_TEST_FILE( input_testNewFrozenVectorDesign )
COPY_TEST_FILE( input_testSolversForUserMatrix-ML )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-2 )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-3 )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-bar )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearThermalOperator-2_HALDEN_clad )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearFickOperator-bar )
COPY_TEST_FILE( input_testTrilinosMueLuSolver-LinearElasticityOperator-1 )
COPY_TEST_FILE( input_testSolversForUserMatrix-MueLu )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-3 )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearThermalOperator-bar )
COPY_TEST_FILE( input_LinearThermalOperator-2-cylinder )
COPY_TEST_FILE( input_LinearThermalOperator-2-shell )
COPY_TEST_FILE( input_LinearThermalOperator-2_HALDEN )
COPY_TEST_FILE( input_LinearThermalOperator-2_HALDEN_clad )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearFickOperator-bar )
COPY_TEST_FILE( input_testBoomerAMGSolver-LinearElasticityOperator-1 )
COPY_TEST_FILE( input_testSolversForUserMatrix-BoomerAMG )

COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-SASolver-HybridGS )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-SASolver-HybridGS-FCG )

COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-IPCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-CG-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-FGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-GCR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-GMRESR-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-TFQMR )

IF ( USE_TRILINOS )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-IPCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-FGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-PetscCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-PetscFGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-ML-PetscBiCGSTAB )

COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-IPCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-FGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-PetscCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-PetscFGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-MueLu-PetscBiCGSTAB )
ENDIF()

IF ( USE_HYPRE )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-IPCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-CG-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-CylMesh-BoomerAMG-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-FGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GCR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-GMRESR-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-TFQMR )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-HypreCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-HypreGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-HypreBiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalPC-HypreBiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-HypreBiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-SASolver-BoomerAMG )

COPY_TEST_FILE( input_testLinearSolvers-PL-AmpMesh-BoomerAMG )
COPY_TEST_FILE( input_testLinearSolvers-PL-CylMesh-BoomerAMG )

IF( USE_PETSC )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscFGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BoomerAMG-PetscBiCGSTAB )
ENDIF()
ENDIF()

IF( USE_PETSC )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-PetscCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-PetscFGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-PetscBiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscFGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-DiagonalSolver-PetscBiCGSTAB )
ENDIF()

COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-FCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-IPCG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-CylMesh-CG )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-GMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-FGMRES )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-BiCGSTAB )
COPY_TEST_FILE( input_testLinearSolvers-LinearThermalRobin-TFQMR )

COPY_TEST_FILE( input_testLinearSolvers-PL-AmpMesh-CG )
COPY_TEST_FILE( input_testLinearSolvers-PL-CylMesh-CG )

COPY_TEST_FILE( input_testPetscSNESSolver-NonlinearThermalOxygenDiffusion-2 )

ADD_AMP_TEST( testSolversWithRhsOptions )

ADD_AMP_TEST( testLinearSolvers-PseudoLaplacian PROCS 1 GPU )
ADD_AMP_TEST( testLinearSolvers-PseudoLaplacian PROCS 2 GPU )
ADD_AMP_TEST( testLinearSolvers-PseudoLaplacian PROCS 4 GPU )

IF ( USE_LIBMESH )
    ADD_AMP_TEST( testLinearSolvers-LinearThermalRobin PROCS 1 GPU )
    ADD_AMP_TEST( testLinearSolvers-LinearThermalRobin PROCS 2 GPU )
    ADD_AMP_TEST( testLinearSolvers-LinearThermalRobin PROCS 4 GPU )

    ADD_AMP_TEST( testSolversForUserMatrix PROCS 1 )
    ADD_AMP_TEST( testSolversForUserMatrix PROCS 2 )
    ADD_AMP_TEST( testSolversForUserMatrix PROCS 4 )

    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-2 PROCS 1 GPU )
    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-2 PROCS 2 GPU )
    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-2 PROCS 4 GPU )

    ADD_AMP_TEST( testLinearSolvers-LinearElasticityOperator-1 PROCS 1 )	
    ADD_AMP_TEST( testLinearSolvers-LinearElasticityOperator-1 PROCS 2 )	
    ADD_AMP_TEST( testLinearSolvers-LinearElasticityOperator-1 PROCS 4 )	

    ADD_AMP_TEST( testLinearSolvers-LinearFickOperator-bar PROCS 1 )	
    ADD_AMP_TEST( testLinearSolvers-LinearFickOperator-bar PROCS 2 )	
    ADD_AMP_TEST( testLinearSolvers-LinearFickOperator-bar PROCS 4 )	

    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-bar PROCS 1 )	
    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-bar PROCS 2 )	
    ADD_AMP_TEST( testLinearSolvers-LinearThermalOperator-bar PROCS 4 )	

    ADD_AMP_TEST( testNonlinearSolvers-NonlinearThermal-1 PROCS 1  )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearThermal-1 PROCS 2 )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearThermal-1 PROCS 4 )

    ADD_AMP_TEST( testNonlinearSolvers-NonlinearMechanics-1 PROCS 1 )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearMechanics-1 PROCS 2 )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearMechanics-1 PROCS 4 )

    ADD_AMP_TEST( testNonlinearSolvers-NonlinearFick PROCS 1 )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearFick PROCS 2 )
    ADD_AMP_TEST( testNonlinearSolvers-NonlinearFick PROCS 4 )

    ADD_AMP_TEST( testLinearSolversWithIdentity GPU )
    ADD_AMP_TEST( testLinearSolversWithIdentity PROCS 2 GPU )
    ADD_AMP_TEST( testLinearSolversWithIdentity PROCS 4 GPU )

    ADD_AMP_PROVISIONAL_TEST( testNonlinearSolvers-NonlinearThermalOxygenDiffusion )

IF( USE_TRILINOS )
    ADD_AMP_TEST( testAMGuniform )
    ADD_AMP_TEST( testMatrixFreeML )
    ADD_AMP_TEST_1_2_4( testLinearRobinBoundaryOperator-1 )
ENDIF()

IF( USE_PETSC AND USE_TRILINOS )
    ADD_AMP_PROVISIONAL_TEST( testMeshRefinementDiffusion-1 )
    ADD_AMP_TEST_1_2_4( testColumnSolver-LinearThermalAndFick-1 )
ENDIF()

ENDIF()

IF( USE_TRILINOS OR USE_PETSC )
    ADD_AMP_TEST( testNonlinearSolversWithIdentity )
ENDIF()

ADD_AMP_TEST( testNewFrozenVectorDesign )

# Add the tests
INITIALIZE_TESTBUILDER()
IF ( USE_LIBMESH AND USE_PETSC AND USE_HYPRE )
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_1 )
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_1 PROCS 2)
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_1 PROCS 4 )
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_2 )
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_2 PROCS 2)
    ADD_TB_TEST( testNonlinearSolvers_NonlinearFickSoret_2 PROCS 4 )
ENDIF()
FINALIZE_TESTBUILDER()

ADD_SUBDIRECTORY( flow )
ADD_SUBDIRECTORY( mechanics )
ADD_SUBDIRECTORY( subchannel )

