# Add the AMP folders and libraries
ADD_AMP_LIBRARY(utils)
ADD_AMP_LIBRARY( IO )
ADD_AMP_LIBRARY( mesh)
ADD_AMP_LIBRARY( geometry)
ADD_AMP_LIBRARY(discretization)
ADD_AMP_LIBRARY(vectors)
ADD_AMP_LIBRARY(matrices)
ADD_AMP_LIBRARY(materials)
ADD_AMP_LIBRARY(operators)
ADD_AMP_LIBRARY(solvers)
ADD_AMP_LIBRARY(time_integrators)
BEGIN_PACKAGE_CONFIG( nek )
IF ( USE_NEK )
    ADD_EXTERNAL_PACKAGE_SUBDIRECTORY( nek ${NEK_EXTERNAL} )
    SET(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -x f77-cpp-input -fcray-pointer" )
    ADD_EXTERNAL_PACKAGE_SUBDIRECTORY( nek/fortran ${NEK_EXTERNAL}/fortran )
    INSTALL_AMP_TARGET ( nek )
ENDIF()
INSTALL_PROJ_LIB()






