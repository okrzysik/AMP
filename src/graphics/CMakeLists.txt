SET( PROJ AMP )
INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )
INCLUDE( "${AMP_SOURCE_DIR}/cmake/configureAMP.cmake" )
SET_AMP_PACKAGE_FLAGS()

BEGIN_PACKAGE_CONFIG( graphics )
INSTALL_AMP_TARGET( graphics )

ADD_SUBDIRECTORY( test )
