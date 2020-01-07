# Add thread pool tests
FILE( COPY "${CMAKE_CURRENT_SOURCE_DIR}/data/commands.txt" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" )
ADD_AMP_TEST( test_atomic )
ADD_AMP_TEST( test_atomic_list )
ADD_AMP_TEST( run_system_tpool 4 commands.txt )
SET_TESTS_PROPERTIES( run_system_tpool--4--commands.txt PROPERTIES PASS_REGULAR_EXPRESSION "echo \"start\"\nstart\n\nsleep 2\; echo \"stop\"\nstop" )
SET_TESTS_PROPERTIES( test_atomic PROPERTIES  FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}"  PROCESSORS 64  RUN_SERIAL 1  COST -500 )
ADD_AMP_TEST_THREAD_MPI( test_thread_pool 1 4 )
ADD_AMP_TEST_THREAD_MPI( test_thread_pool 2 4 )
ADD_AMP_TEST_THREAD_MPI( test_thread_pool 4 4 )
SET_TESTS_PROPERTIES( test_thread_pool_4threads PROPERTIES  FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}"  PROCESSORS 4  RUN_SERIAL 1  COST -100 )
IF ( USE_MPI )
    IF ( 8 LESS_EQUAL ${TEST_MAX_PROCS} )
        SET_TESTS_PROPERTIES( test_thread_pool_2procs_4threads PROPERTIES  FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}"  PROCESSORS 8  RUN_SERIAL 1  COST -200 )
    ENDIF()
    IF ( 16 LESS_EQUAL ${TEST_MAX_PROCS} )
        SET_TESTS_PROPERTIES( test_thread_pool_4procs_4threads PROPERTIES  FAIL_REGULAR_EXPRESSION "${TEST_FAIL_REGULAR_EXPRESSION}"  PROCESSORS 16  RUN_SERIAL 1  COST -400 )
    ENDIF()
ENDIF()

# Add affinity tests
COPY_TEST_FILE( input_print_affinity_1 )
ADD_AMP_TEST( print_thread_affinities input_print_affinity_1 )
ADD_AMP_PROVISIONAL_TEST( test_async_communication )


