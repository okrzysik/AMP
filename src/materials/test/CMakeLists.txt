############################################################################
#  Configure tests
############################################################################

ADD_AMP_TEST ( testMaterial )
ADD_AMP_TEST ( test_UO2_MSRZC_09 )
ADD_AMP_TEST ( testWaterLibrary )
ADD_AMP_TEST ( test_Independent )
ADD_AMP_TEST ( examiner )
ADD_AMP_TEST ( testExtraParams2Eval )
ADD_AMP_TEST ( test_FixedClad )
ADD_AMP_TEST ( test_FixedFuel )

COPY_TEST_FILE( inputExaminer )
