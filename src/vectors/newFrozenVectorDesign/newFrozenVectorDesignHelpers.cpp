
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include "AMP/vectors/MultiVector.h"

namespace AMP::LinearAlgebra {

AMP::LinearAlgebra::Vector::shared_ptr
subsetExceptForVariable( AMP::LinearAlgebra::Vector::shared_ptr inVec,
                         std::shared_ptr<AMP::LinearAlgebra::Variable> var )
{
    auto outVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
        AMP::LinearAlgebra::MultiVector::create( "MultiVariable", AMP_COMM_SELF ) );
    auto castedInVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( inVec );
    auto vec2skip    = inVec->subsetVectorForVariable( var );
    for ( auto vec : *castedInVec )
        if ( vec != vec2skip ) {
            outVec->addVector( vec );
        }

    return outVec;
}

AMP::LinearAlgebra::Vector::shared_ptr joinVectors( AMP::LinearAlgebra::Vector::shared_ptr vec1,
                                                    AMP::LinearAlgebra::Vector::shared_ptr vec2 )
{
    auto outVec = std::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
        AMP::LinearAlgebra::MultiVector::create( "MultiVariable", AMP_COMM_SELF ) );
    outVec->addVector( vec1 );
    outVec->addVector( vec2 );
    return outVec;
}
} // namespace AMP::LinearAlgebra
