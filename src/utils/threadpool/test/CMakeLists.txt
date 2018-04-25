# Add thread pool tests
ADD_AMP_TEST( test_atomic )
ADD_AMP_TEST( test_atomic_list )
SET_TESTS_PROPERTIES( test_atomic PROPERTIES  FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}"  PROCESSORS 64  RUN_SERIAL 1  COST -500 )
ADD_AMP_TEST_THREAD_MPI( test_thread_pool 1 4 )
IF ( USE_MPI )
    ADD_AMP_TEST_THREAD_MPI( test_thread_pool 2 4 )
    ADD_AMP_TEST_THREAD_MPI( test_thread_pool 4 4 )
ENDIF()

# Add affinity tests
COPY_TEST_FILE( input_print_affinity_1 )
ADD_AMP_TEST( print_thread_affinities input_print_affinity_1 )
ADD_AMP_PROVISIONAL_TEST( test_async_communication )
