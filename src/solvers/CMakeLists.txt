INCLUDE (macros)

BEGIN_PACKAGE_CONFIG ( solvers )
ADD_PACKAGE_SUBDIRECTORY ( newFrozenVectorDesign )
IF( USES_DENDRO )
    ADD_PACKAGE_SUBDIRECTORY ( contact )
ENDIF()
INSTALL_AMP_TARGET ( solvers )

ADD_SUBDIRECTORY ( test )
ADD_SUBDIRECTORY ( newFrozenVectorDesign/test )

