#include "vectors/trilinos/NativeThyraVector.h"

#include "Thyra_SpmdVectorBase_def.hpp"
#include "Thyra_DefaultSpmdVector_def.hpp"


namespace AMP {
namespace LinearAlgebra {

NativeThyraVector::NativeThyraVector ( VectorParameters::shared_ptr in_params ): 
    NativeVector(), 
    ThyraVector(), 
    VectorEngine()
{ 
    boost::shared_ptr<NativeThyraVectorParameters> params = 
        boost::dynamic_pointer_cast<NativeThyraVectorParameters>(in_params);
    AMP_ASSERT(params!=NULL);
    AMP_ASSERT(!params->d_comm.isNull());
    AMP_ASSERT(params->d_InVec.get()!=NULL);
    Thyra::Ordinal dim = params->d_InVec->space()->dim();
    AMP_ASSERT(params->d_comm.sumReduce(params->d_local)==static_cast<size_t>(dim));
    boost::shared_ptr<CommunicationListParameters> communicationListParams( new CommunicationListParameters() );
    communicationListParams->d_comm = params->d_comm;
    communicationListParams->d_localsize = params->d_local;
    d_CommList = boost::shared_ptr<CommunicationList>( new CommunicationList(communicationListParams ) );
    d_DOFManager = boost::shared_ptr<Discretization::DOFManager>( new Discretization::DOFManager(params->d_local,params->d_comm) );
    d_local = params->d_local;
    d_thyraVec = params->d_InVec;
    d_pVariable = params->d_var;
}


NativeThyraVector::~NativeThyraVector ()
{
}


Vector::shared_ptr NativeThyraVector::cloneVector(const Variable::shared_ptr var ) const 
{ 
    boost::shared_ptr<NativeThyraVectorParameters> params( new NativeThyraVectorParameters() );
    params->d_InVec = d_thyraVec->clone_v();
    params->d_local = d_local;
    params->d_comm = getComm();
    params->d_var = d_pVariable;
    return boost::shared_ptr<NativeThyraVector>( new NativeThyraVector(params) );
}


void NativeThyraVector::putRawData ( double *in )
{
    AMP_ERROR( "not implemented" );
}


void* NativeThyraVector::getRawDataBlockAsVoid ( size_t i )
{ 
    Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    Thyra::DefaultSpmdVector<double>* spmdVector = 
        dynamic_cast<Thyra::DefaultSpmdVector<double>*>(ptr);
    if ( spmdVector!=NULL ) {
        if ( i!=0 )
            AMP_ERROR("Invalid block");
        return spmdVector->getPtr();
    }
    AMP_ERROR( "not finished" );
    return NULL;
}


const void* NativeThyraVector::getRawDataBlockAsVoid ( size_t i ) const
{ 
    const Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    const Thyra::DefaultSpmdVector<double>* spmdVector = 
        dynamic_cast<const Thyra::DefaultSpmdVector<double>*>(ptr);
    if ( spmdVector!=NULL ) {
        if ( i!=0 )
            AMP_ERROR("Invalid block");
        return spmdVector->getPtr();
    }
    return NULL;
}


size_t NativeThyraVector::numberOfDataBlocks () const 
{ 
    const Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    const Thyra::DefaultSpmdVector<double>* spmdVector = 
        dynamic_cast<const Thyra::DefaultSpmdVector<double>*>(ptr);
    if ( spmdVector!=NULL ) {
        return 1;
    }
    AMP_ERROR( "not finished" );
    return 1;
}


size_t NativeThyraVector::sizeOfDataBlock ( size_t i ) const
{
    const Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    const Thyra::DefaultSpmdVector<double>* spmdVector = 
        dynamic_cast<const Thyra::DefaultSpmdVector<double>*>(ptr);
    if ( spmdVector!=NULL ) {
        return d_local;
    }
    AMP_ERROR( "not finished" );
    return d_local;
}


}
}

