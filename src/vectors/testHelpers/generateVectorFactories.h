// Generate the desired vector factory
#ifndef included_test_generateVectorFactory
#define included_test_generateVectorFactory


#include "vectors/testHelpers/VectorFactory.h"
#include "vectors/testHelpers/VectorTests.h"


namespace AMP {
namespace LinearAlgebra {


/** \brief Generate the desired vector factory
  * \details  This will generate the desired vector factory
  * \param[in] factory      The vector factory to generate:
        SimpleVectorFactory<15,false,double>
        SimpleVectorFactory<45, true,double>
        SimpleVectorFactory<45,true,double,openmp,gpu>
        ArrayVectorFactory<4,10,false,double>
        ArrayVectorFactory<4,10,true,double>
        SimplePetscNativeFactory
        SimpleManagedVectorFactory<ManagedEpetraVector>
        MultiVectorFactory<SMEVFactory, 1, SNPVFactory, 1> MVFactory1;
        MultiVectorFactory<SMEVFactory, 3, SNPVFactory, 2> MVFactory2;
        MultiVectorFactory<MVFactory1, 2, MVFactory2, 2> MVFactory3;
        MultiVectorFactory<SimpleVectorFactory<15,false>,1,SMEVFactory,1> MVFactory1;
        MultiVectorFactory<SimpleVectorFactory<15,false>,3,SMEVFactory,2> MVFactory2;
        MultiVectorFactory<MVFactory1, 2, MVFactory2, 2> MVFactory3;
  *
  */
AMP::shared_ptr<VectorFactory> generateVectorFactory( const std::string &factory );
}
}

/// \endcond

#endif
