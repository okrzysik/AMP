cmake                                           ^
      -G "NMake Makefiles"                      ^
      -D AMP_DATA:PATH=T:/AMP/AMP-Data          ^
      -D COMPILE_MODE:STRING=debug              ^
      -D USE_FORTRAN=0                          ^
      -D USE_EXT_BOOST=1                        ^
         -D BOOST_DIRECTORY=T:/AMP/boost_1_53_0 ^
      -D USE_EXT_BLAS=0                         ^
      -D USE_EXT_LAPACK=0                       ^
      -D USE_EXT_TRILINOS=0                     ^
      -D USE_EXT_PETSC=0                        ^
      -D USE_EXT_LIBMESH=0                      ^
      -D USE_EXT_SUNDIALS=0                     ^
      -D USE_EXT_SILO=0                         ^
      -D USE_EXT_HDF5=0                         ^
      -D USE_EXT_HYPRE=0                        ^
      -D USE_EXT_X11=0                          ^
      -D USE_EXT_MPI=0                          ^
      -D USE_AMP_UTILS=1                        ^
      -D USE_AMP_MESH=1                         ^
      -D USE_AMP_DISCRETIZATION=1               ^
      -D USE_AMP_VECTORS=1                      ^
      -D USE_AMP_MATRICES=1                     ^
      -D USE_AMP_MATERIALS=1                    ^
      -D USE_AMP_OPERATORS=0                    ^
      -D USE_AMP_SOLVERS=0                      ^
      -D USE_AMP_TIME_INTEGRATORS=0             ^
      ../../AMP





