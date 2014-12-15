INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE( "${AMP_SOURCE_DIR}/cmake/configureAMP.cmake" )
SET_AMP_PACKAGE_FLAGS()

BEGIN_PACKAGE_CONFIG ( utils )
INSTALL_AMP_TARGET ( utils )

ADD_SUBDIRECTORY ( test )
