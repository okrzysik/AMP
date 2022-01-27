BEGIN_EXAMPLE_CONFIG( utils )

# Copy meshes and input files
COPY_EXAMPLE_DATA_FILE( input_readInputFile )

# Add the examples:
#    ADD_AMP_EXAMPLE( exampleName  N_procs  inputArgs )
# Note: each test executable must be added at least once to ensure it compiles

ADD_AMP_EXAMPLE( ampHelloWorld 1 )

IF ( USE_MPI )

ADD_AMP_EXAMPLE( ampHelloWorld 2 )
ADD_AMP_EXAMPLE( ampHelloWorld 4 )

ADD_AMP_EXAMPLE( ampLogging 2 )

ENDIF()

ADD_AMP_EXAMPLE( readInputFile 1 input_readInputFile)


# Install the examples
INSTALL_AMP_EXAMPLE( utils )
