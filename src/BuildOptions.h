//---------------------------------------------------------------------------//
// BuildOptions.h
// Include file for doxygen providing the build options
//---------------------------------------------------------------------------//


/*! \page BuildOptions  AMP configure options
  This page provides a list of the different compile options for AMP.  AMP is configured using cmake
  and the order of the options for cmake does not matter.

  \par AMP internal package options
  \ref AMP::Utilities "UTILS": The UTILS package provides all of the utility
      functions used throughout AMP and is required to build AMP. <BR>
  \ref AMP::Mesh "MESH": The Mesh package provides all of the mesh interface routines and is
      required for a number of additional packages. <BR>
  \ref AMP::Discretization "DISCRETIZATION": The DISCRETIZATION package provides classes for
  managing DOFs
      (Degrees of Freedom) and discretization.  It requires MESH. <BR>
  \ref AMP::LinearAlgebra "VECTORS": The VECTORS package provides the vector classes
      in LinearAlgebra that are used within AMP<BR>
  \ref AMP::LinearAlgebra "MATRICES": The MATRICES package provides the matrix classes
      in LinearAlgebra that are used within AMP.  It requires VECTORS. <BR>
  \ref AMP::Materials "MATERIALS": The MATERIALS package provides the materials interface and
  requires VECTORS. <BR>
  \ref AMP::Operator "OPERATORS": The OPERATORS package provides the interface and all operators
  within AMP.
      It requires MESH, VECTORS, MATRICES, and the external package LIBMESH. <BR>
  \ref AMP::Solver "SOLVERS": The SOLVERS package provides routines for solving and requires
  OPERATORS. <BR>
  \ref AMP::TimeIntegrator "TIME_INTEGRATORS": The TIME_INTEGRATORS package provides routines for
  time integration and
  requires SOLVERS. <BR>

  \par AMP external package options
  AMP is designed to leverage a number of external packages.
  While most external packages are optional, some functionality will be lost if compiling without
  a given package.  All external packages can be included or disabled using the "-D
  USE_EXT_PACKAGE="
  flags.  Unless otherwise noted, all external packages take a second cmake argument "-D
  PACKAGE_DIRECTORY="
  specifying the directory where the given package is installs.  The current list of external
  packages is: <BR>
  MPI: MPI provides parallel capabilities for AMP.  If used, there are a number of additional
  optional flags
     that are used to configure MPI: <BR>
        &nbsp;&nbsp;&nbsp; "-D USE_MPI_FOR_SERIAL_TESTS=" indicates that we want to use mpi to run
  the serial tests in
  ctest.  Default is false. <BR>
        &nbsp;&nbsp;&nbsp; "-D MPIEXEC_NUMPROC_FLAG:STRING=" the flag to specify the flag mpi uses
  to specify the number
  of processors.  Default is "-np"<BR>
  LAPACK: Blas libraries that are used by many external packages. <BR>
  BLAS: Blas libraries that are used by many external packages. <BR>
  X11: X11 Libraries that may be required by some external packages. <BR>
  TRILINOS: Trilinos is used throughout AMP to provide additional functionality.  It is used but not
     required within LinearAlgebra, Operators and Solvers.  Compiling without Trilinos will result
  in reduced
  functionality.  <BR>
  PETSC: Petsc is used throughout AMP to provide additional functionality.  It is used but not
     required within LinearAlgebra, Operators and Solvers.  Compiling without Petsc will result in
  reduced
  functionality.  <BR>
  LIBMESH: LibMesh is required for the Operators package.  It is currently used within Mesh for
  additional
  functionality.  <BR>
  STKMESH: stk::mesh is currently used within Mesh for additional functionality.  <BR>
  SUNDIALS: <BR>
  HYPRE:<BR>
  NEK:<BR>
  MOAB:<BR>
  DENDRO:<BR>
  SILO:<BR>
  HDF5:<BR>
  NETCDF:<BR>

  \par Additional compile flags
  "-D COMPILE_MODE=":  Required string indicating if we are building in "debug" or "optimized mode".
  <BR>
  "-D AMP_DATA=":  Required path to the amp data directory.  <BR>
  "-D CMAKE_C_COMPILER="  Optional flag indicating the C compiler to use.  <BR>
  "-D CMAKE_CXX_COMPILER="  Optional flag indicating the CXX compiler to use.  <BR>
  "-D CMAKE_Fortran_COMPILER="  Optional flag indicating the Fortran compiler to use.  <BR>
  "-D CFLAGS="  Optional flag indicating additional flags to pass to the C compiler <BR>
  "-D CXXFLAGS="  Optional flag indicating additional flags to pass to the CXX compiler <BR>
  "-D FFLAGS="  Optional flag indicating additional flags to pass to the fortran compiler <BR>
  "-D LDLIBS="  Optional flag indicating additional flags to pass to the linker <BR>
  "-D TEST_MAX_PROCS="  Optional flag indicating the maximum number of processors for the tests
  (default: 32) <BR>
  "-D USE_STATIC="  Optional flag indicating if we are only using static libraries (default: false)
  <BR>
  "-D USE_FORTRAN="  Optional flag indicating if we want to include the fortran compiler (default:
  true) <BR>

*/
