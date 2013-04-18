INCLUDE (${AMP_SOURCE_DIR}/cmake/macros.cmake)
INCLUDE (${AMP_SOURCE_DIR}/cmake/libraries.cmake)
SET_AMP_PACKAGE_FLAGS()

BEGIN_PACKAGE_CONFIG ( time_integrators )
ADD_PACKAGE_SUBDIRECTORY ( oxide )

IF ( USE_EXT_SUNDIALS )
    ADD_PACKAGE_SUBDIRECTORY ( sundials )
ENDIF()

INSTALL_AMP_TARGET ( time_integrators )

ADD_SUBDIRECTORY ( test )

