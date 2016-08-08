BEGIN_PACKAGE_CONFIG ( matrices )
IF ( USE_EXT_TRILINOS )
    ADD_PACKAGE_SUBDIRECTORY ( trilinos )
ENDIF()
IF ( USE_EXT_PETSC AND USE_EXT_TRILINOS )
    ADD_PACKAGE_SUBDIRECTORY ( petsc )
ELSEIF ( USE_EXT_PETSC )
    # Note: currently petsc matricies require trilinos, this needs to be fixed
    MESSAGE( STATUS "PETSc matrix support currently requires Trilinos" )
ENDIF()
INSTALL_AMP_TARGET ( matrices )

ADD_SUBDIRECTORY ( test )
