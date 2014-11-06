#include "vectors/VectorSelector.h"
#include "vectors/MultiVector.h"
#include "vectors/VectorSelector.h"

#include "utils/AMPManager.h"
#include "utils/UnitTest.h"



template <typename T>
class StridedVectorFactory
{
public:

    static AMP::LinearAlgebra::Variable::shared_ptr  getVariable()
    {
        return AMP::LinearAlgebra::Variable::shared_ptr ();
    }

    static AMP::LinearAlgebra::Vector::shared_ptr   getVector()
    {
        AMP::shared_ptr<typename T::vector>  vec = T::getVector();
        AMP::LinearAlgebra::VS_Stride criterion = AMP::LinearAlgebra::VS_Stride(1,3);
        AMP::LinearAlgebra::Vector::shared_ptr  vec_select = vec->select( criterion, "thirds" );
        size_t N1 = vec->getGlobalSize();
        size_t N2 = vec_select->getGlobalSize();
        AMP_ASSERT(N1/3==N2);
        return vec_select;
    }
};


