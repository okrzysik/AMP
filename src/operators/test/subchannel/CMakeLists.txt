INCLUDE( "${AMP_SOURCE_DIR}/cmake/macros.cmake" )

ADD_AMP_TEST(testSubchannelHelpers)
COPY_TEST_DATA_FILE ( input_testSubchannelUtilities)

ADD_AMP_TEST_1_2_4 (testSubchannelPhysicsModel)
COPY_TEST_DATA_FILE ( input_testSubchannelPhysicsModel)

ADD_AMP_TEST(testSubchannelTwoEqNonlinearOperator)
COPY_TEST_DATA_FILE ( input_testSubchannelTwoEqNonlinearOperator)

ADD_AMP_TEST(testSubchannelTwoEqLinearOperator)
COPY_TEST_DATA_FILE ( input_testSubchannelTwoEqLinearOperator)

ADD_AMP_TEST(testSubchannelFourEqNonlinearOperator)
COPY_TEST_DATA_FILE ( input_testSubchannelFourEqNonlinearOperator)

ADD_AMP_TEST(testSubchannelFourEqLinearOperator)
COPY_TEST_DATA_FILE ( input_testSubchannelFourEqLinearOperator)




