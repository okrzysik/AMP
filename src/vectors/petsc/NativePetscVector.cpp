#include "AMP/vectors/petsc/NativePetscVector.h"
#include "AMP/vectors/petsc/NativePetscVectorOperations.h"

#include "AMP/vectors/petsc/ManagedPetscVector.h"

#include "petsc.h"
#include "petsc/private/vecimpl.h"
#include "petscvec.h"

namespace AMP {
namespace LinearAlgebra {

NativePetscVector::NativePetscVector( VectorParameters::shared_ptr in_params )
    : Vector()
{
    d_VectorOps    = std::make_shared<NativePetscVectorOperations>();
    d_VectorData   = new NativePetscVectorData(in_params);
}

NativePetscVector::NativePetscVector( std::shared_ptr<VectorData> data )
  : Vector(), d_VectorDataSP{data}
{
    d_VectorOps    = std::make_shared<NativePetscVectorOperations>();
    d_VectorData   = data.get();
}

NativePetscVector::~NativePetscVector()
{
  if (d_VectorData) delete d_VectorData;
}

Vector::shared_ptr NativePetscVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto data = d_VectorData->cloneData();
    auto retVal = std::make_shared<NativePetscVector>( data );
    retVal->setVariable( var );
    return retVal;
}

std::shared_ptr<ParameterBase> NativePetscVector::getParameters()
{ return dynamic_cast<NativePetscVector*>(d_VectorData)->getParameters(); }

void NativePetscVector::aliasVector( Vector & ) { AMP_ERROR( "not implemented" ); }

void NativePetscVector::swapVectors( Vector &other )
{
  d_VectorData->swapData(*other.getVectorData());
}

} // namespace LinearAlgebra
} // namespace AMP
