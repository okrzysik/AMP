############################################################################
#  Configure tests
############################################################################

# Copy test inputs


# Copy test data


# Add the tests
ADD_AMP_PROVISIONAL_TEST( writeBox )
ADD_AMP_PROVISIONAL_TEST( write2elementMesh )
ADD_AMP_PROVISIONAL_TEST( write3pointConstrainedBox )
ADD_AMP_PROVISIONAL_TEST( writePlateWithHoleMesh )
ADD_AMP_TEST( write7elementMesh-1 out7elementMesh-1 )
ADD_AMP_TEST( write7elementMesh-2 out7elementMesh-2 )
ADD_AMP_TEST( writeDistortedElementMesh outDistortedElementMesh )
ADD_AMP_TEST( writeDispValsForPatchTest-2 outDispValsForPatchTest-2 )




