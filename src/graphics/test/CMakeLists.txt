# Graphics tests

IF ( USE_EXT_OGRE )
    COPY_TEST_FILE( resources.cfg )
    COPY_TEST_FILE( plugins.cfg )
    ADD_AMP_PROVISIONAL_TEST( test_ogre_tutorial1 )
ENDIF()

