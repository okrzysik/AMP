# Graphics tests
IF ( USE_EXT_OGRE )
    CONFIGURE_FILE( plugins.cfg "${CMAKE_CURRENT_BINARY_DIR}/${plugins.cfg}" )
    CONFIGURE_FILE( resources.cfg "${CMAKE_CURRENT_BINARY_DIR}/${resources.cfg}" )
    ADD_AMP_PROVISIONAL_TEST( test_ogre_tutorial1 )
ENDIF()

