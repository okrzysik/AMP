BEGIN_EXAMPLE_CONFIG( mesh )

# Copy meshes and input files
COPY_EXAMPLE_DATA_FILE( input_createMesh-cube )
COPY_EXAMPLE_DATA_FILE( input_createMesh-2d )
COPY_EXAMPLE_DATA_FILE( input_createMesh-3d )
COPY_EXAMPLE_DATA_FILE( input_createMesh-Multimesh )

# Add the examples:
#    ADD_AMP_EXAMPLE( exampleName  N_procs  inputArgs )
# Note: each test executable must be added at least once to ensure it compiles
ADD_AMP_EXAMPLE( createMesh  1  input_createMesh-cube )
ADD_AMP_EXAMPLE( createMesh  1  input_createMesh-2d )
ADD_AMP_EXAMPLE( createMesh  1  input_createMesh-3d )
ADD_AMP_EXAMPLE( createMesh  1  input_createMesh-Multimesh )


# Install the examples
INSTALL_AMP_EXAMPLE( mesh )

