
#include "AMP/vectors/newFrozenVectorDesign/newFrozenVectorDesignHelpers.h"

#include "AMP/vectors/MultiVector.h"

namespace AMP {
namespace LinearAlgebra {

AMP::LinearAlgebra::Vector::shared_ptr
subsetExceptForVariable( AMP::LinearAlgebra::Vector::shared_ptr inVec,
                         AMP::LinearAlgebra::Variable::shared_ptr var )
{
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> outVec =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
            AMP::LinearAlgebra::MultiVector::create( "MultiVariable", AMP_COMM_SELF ) );
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> castedInVec =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>( inVec );
    AMP::LinearAlgebra::Vector::shared_ptr vec2skip = inVec->subsetVectorForVariable( var );
    auto curr                                       = castedInVec->beginVector();
    auto end                                        = castedInVec->endVector();
    for ( ; curr != end; curr++ ) {
        if ( ( *curr ) != vec2skip ) {
            outVec->addVector( *curr );
        }
    }

    return outVec;
}

AMP::LinearAlgebra::Vector::shared_ptr joinVectors( AMP::LinearAlgebra::Vector::shared_ptr vec1,
                                                    AMP::LinearAlgebra::Vector::shared_ptr vec2 )
{
    AMP::shared_ptr<AMP::LinearAlgebra::MultiVector> outVec =
        AMP::dynamic_pointer_cast<AMP::LinearAlgebra::MultiVector>(
            AMP::LinearAlgebra::MultiVector::create( "MultiVariable", AMP_COMM_SELF ) );
    outVec->addVector( vec1 );
    outVec->addVector( vec2 );
    return outVec;
}
} // namespace LinearAlgebra
} // namespace AMP
