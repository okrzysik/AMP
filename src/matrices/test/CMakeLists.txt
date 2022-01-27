############################################################################
#  Configure tests
############################################################################

# Function to add a serial or 1,2,4 processor test dependent on the existence of a parallel matrix
FUNCTION( ADD_AMP_TEST_124_MATRIX EXEFILE ${ARGN} )
    IF ( USE_TRILINOS )
        # We have a parallel matrix
        ADD_AMP_TEST_1_2_4( ${EXEFILE} ${ARGN} )
    ELSE()
        # We do not have a parallel matrix
        ADD_AMP_TEST( ${EXEFILE} ${ARGN} )
    ENDIF()
ENDFUNCTION()


# Copy meshes
IF ( AMP_DATA )
    COPY_MESH_FILE( pellet_1x.e )
    COPY_MESH_FILE( pellet_lo_res.e )
    COPY_MESH_FILE( clad_1x_1pellet.e )
    COPY_MESH_FILE( clad_1x_2pellet.e )
    COPY_MESH_FILE( cube343H27.e )
    COPY_MESH_FILE( mpcMesh-1 )
ENDIF()


# Add matrix test
ADD_AMP_TEST_1_2_4( test_Matrix )
ADD_AMP_TEST_1_2_4( test_RectangularMatricesMultiplication )

IF( USE_LIBMESH AND AMP_DATA )
  ADD_AMP_TEST_124_MATRIX( testMatCreate )
  COPY_TEST_FILE( input_testMatCreate-1 )
  ADD_AMP_TEST_1_2_4( testMatMultiply )
ENDIF()

ADD_SUBDIRECTORY ( extra_tests )

