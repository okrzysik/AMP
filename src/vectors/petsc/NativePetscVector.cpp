#include "AMP/vectors/petsc/NativePetscVector.h"
#include "AMP/vectors/petsc/NativePetscVectorOperations.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector( Vec v, bool deleteable, AMP_MPI comm ) : Vector()
{
    d_VectorOps  = std::make_shared<NativePetscVectorOperations>();
    d_VectorData = std::make_shared<NativePetscVectorData>( v, deleteable, comm );
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( d_VectorData->getLocalSize(),
                                                                      d_VectorData->getComm() );
}

NativePetscVector::NativePetscVector( std::shared_ptr<VectorData> data ) : Vector()
{
    d_VectorData = data;
    d_VectorOps  = std::make_shared<NativePetscVectorOperations>();
    d_DOFManager =
        std::make_shared<AMP::Discretization::DOFManager>( data->getLocalSize(), data->getComm() );
}

NativePetscVector::~NativePetscVector() {}

Vector::shared_ptr NativePetscVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto data   = d_VectorData->cloneData();
    auto retVal = std::make_shared<NativePetscVector>( data );
    retVal->setVariable( var );
    return retVal;
}

void NativePetscVector::swapVectors( Vector &other )
{
    d_VectorData->swapData( *other.getVectorData() );
}

} // namespace LinearAlgebra
} // namespace AMP
