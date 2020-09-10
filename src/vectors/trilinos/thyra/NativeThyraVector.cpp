#include "AMP/vectors/trilinos/thyra/NativeThyraVector.h"
#include "AMP/vectors/trilinos/thyra/NativeThyraVectorData.h"
#include "AMP/vectors/trilinos/thyra/NativeThyraVectorOperations.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


// Trilinos includes
DISABLE_WARNINGS
//#include "Thyra_SpmdVectorBase_def.hpp"
#include "Thyra_DefaultSpmdVector_def.hpp"
#include "Thyra_VectorStdOps_def.hpp"
#include "Trilinos_version.h"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {


/************************************************************************
 * Constructors                                                          *
 ************************************************************************/
NativeThyraVector::NativeThyraVector( VectorParameters::shared_ptr in_params ) : Vector()
{
    d_VectorOps = std::make_shared<NativeThyraVectorOperations>();

    setVectorData( std::make_shared<NativeThyraVectorData>( in_params ) );

    auto params = std::dynamic_pointer_cast<NativeThyraVectorParameters>( in_params );
    AMP_ASSERT( params != nullptr );
    AMP_ASSERT( !params->d_comm.isNull() );

    d_DOFManager = std::make_shared<Discretization::DOFManager>( params->d_local, params->d_comm );
    d_pVariable  = params->d_var;
}

NativeThyraVector::NativeThyraVector( std::shared_ptr<VectorData> data ) : Vector()
{
    setVectorData( data );
    d_VectorOps = std::make_shared<NativeThyraVectorOperations>();
    d_DOFManager =
        std::make_shared<AMP::Discretization::DOFManager>( data->getLocalSize(), data->getComm() );
}

/************************************************************************
 * Destructor                                                            *
 ************************************************************************/
NativeThyraVector::~NativeThyraVector() {}

/************************************************************************
 * Vector functions                                                      *
 ************************************************************************/
Vector::shared_ptr NativeThyraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto data   = d_VectorData->cloneData();
    auto retVal = std::make_shared<NativeThyraVector>( data );
    retVal->setVariable( var );
    return retVal;
}

void NativeThyraVector::swapVectors( Vector & ) { AMP_ERROR( "not implemented" ); }

} // namespace LinearAlgebra
} // namespace AMP
