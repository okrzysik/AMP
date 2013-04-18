INCLUDE (${AMP_SOURCE_DIR}/cmake/macros.cmake)
INCLUDE (${AMP_SOURCE_DIR}/cmake/libraries.cmake)
SET_AMP_PACKAGE_FLAGS()

BEGIN_PACKAGE_CONFIG ( materials )
INSTALL_AMP_TARGET ( materials )

ADD_SUBDIRECTORY ( test )
