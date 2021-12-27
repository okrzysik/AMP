############################################################################
#  Configure tests
############################################################################

# Copy inputs
COPY_TEST_FILE( input_geometry_zoo )

# Add tests
ADD_AMP_PROVISIONAL_TEST( test_Geometry )
#ADD_AMP_TEST( test_Geometry input_geometry_zoo )
ADD_AMP_TEST( test_GeometryHelpers )

