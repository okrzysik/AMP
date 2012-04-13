INCLUDE (macros)



# Add the AMP folders and libraries
IF ( USE_AMP_UTILS )
    ADD_AMP_LIBRARY (utils)
ENDIF()
IF ( USE_AMP_MESH )
    ADD_AMP_LIBRARY (ampmesh)
ENDIF()
IF ( USE_AMP_DISCRETIZATION )
    ADD_AMP_LIBRARY (discretization)
ENDIF()
IF ( USE_AMP_VECTORS )
    ADD_AMP_LIBRARY (vectors)
ENDIF()
IF ( USE_AMP_MATRICIES )
    ADD_AMP_LIBRARY (matrices)
ENDIF()
IF ( USE_AMP_MATERIALS )
    ADD_AMP_LIBRARY (materials)
ENDIF()
IF ( USE_AMP_OPERATORS )
    ADD_AMP_LIBRARY (operators)
ENDIF()
IF ( USE_AMP_SOLVERS )
    ADD_AMP_LIBRARY (solvers)
ENDIF()
IF ( USE_AMP_TIME_INTEGRATORS )
    ADD_AMP_LIBRARY (time_integrators)
ENDIF()
IF ( USE_NEK )
    BEGIN_PACKAGE_CONFIG( nek )
    ADD_EXTERNAL_PACKAGE_SUBDIRECTORY( nek ${NEK_EXTERNAL} )
    SET(CMAKE_Fortran_FLAGS " ${CMAKE_Fortran_FLAGS} -x f77-cpp-input -fcray-pointer" )
    ADD_EXTERNAL_PACKAGE_SUBDIRECTORY( nek/fortran ${NEK_EXTERNAL}/fortran )
    INSTALL_AMP_TARGET ( nek )
ENDIF()







