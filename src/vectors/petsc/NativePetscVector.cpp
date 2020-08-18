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
    d_VectorDataSP = std::make_shared<NativePetscVectorData>(in_params);
    d_VectorData   = d_VectorDataSP.get();
    
    auto npvParams = std::dynamic_pointer_cast<NativePetscVectorParameters>( in_params );   
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( npvParams->d_localsize,
                                                                      npvParams->d_Comm );
}

NativePetscVector::NativePetscVector( std::shared_ptr<VectorData> data )
  : Vector(), d_VectorDataSP{data}
{
    d_VectorOps    = std::make_shared<NativePetscVectorOperations>();
    d_VectorData   = data.get();
    d_DOFManager = std::make_shared<AMP::Discretization::DOFManager>( data->getLocalSize(),
                                                                      data->getComm() );
}

NativePetscVector::~NativePetscVector()
{
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
void NativePetscVector::assemble( )
{
  dynamic_cast<NativePetscVectorData*>(d_VectorData)->assemble();
}

} // namespace LinearAlgebra
} // namespace AMP
