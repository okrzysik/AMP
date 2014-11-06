#include "vectors/trilinos/thyra/NativeThyraVector.h"
#include "vectors/trilinos/thyra/ThyraVectorWrapper.h"

//#include "Thyra_SpmdVectorBase_def.hpp"
#include "Thyra_DefaultSpmdVector_def.hpp"


namespace AMP {
namespace LinearAlgebra {


/************************************************************************
* Constructors                                                          *
************************************************************************/
NativeThyraVector::NativeThyraVector ( VectorParameters::shared_ptr in_params ): 
    NativeVector(), 
    ThyraVector(), 
    VectorEngine()
{ 
    AMP::shared_ptr<NativeThyraVectorParameters> params = 
        AMP::dynamic_pointer_cast<NativeThyraVectorParameters>(in_params);
    AMP_ASSERT(params!=NULL);
    AMP_ASSERT(!params->d_comm.isNull());
    AMP_ASSERT(params->d_InVec.get()!=NULL);
    Thyra::Ordinal dim = params->d_InVec->space()->dim();
    AMP_ASSERT(params->d_comm.sumReduce(params->d_local)==static_cast<size_t>(dim));
    AMP::shared_ptr<CommunicationListParameters> communicationListParams( new CommunicationListParameters() );
    communicationListParams->d_comm = params->d_comm;
    communicationListParams->d_localsize = params->d_local;
    d_CommList = AMP::shared_ptr<CommunicationList>( new CommunicationList(communicationListParams ) );
    d_DOFManager = AMP::shared_ptr<Discretization::DOFManager>( new Discretization::DOFManager(params->d_local,params->d_comm) );
    d_local = params->d_local;
    d_thyraVec = params->d_InVec;
    d_pVariable = params->d_var;
}


/************************************************************************
* Destructor                                                            *
************************************************************************/
NativeThyraVector::~NativeThyraVector ()
{
}


/************************************************************************
* Vector functions                                                      *
************************************************************************/
Vector::shared_ptr NativeThyraVector::cloneVector(const Variable::shared_ptr var ) const 
{ 
    AMP::shared_ptr<NativeThyraVectorParameters> params( new NativeThyraVectorParameters() );
    params->d_InVec = d_thyraVec->clone_v();
    params->d_local = d_local;
    params->d_comm = getComm();
    params->d_var = var;
    return AMP::shared_ptr<NativeThyraVector>( new NativeThyraVector(params) );
}


void NativeThyraVector::putRawData( const double *in )
{
    size_t i=0;
    for (size_t b=0; b<numberOfDataBlocks(); b++) {
        double *data = reinterpret_cast<double*>(getRawDataBlockAsVoid(b));
        for (size_t j=0; j<sizeOfDataBlock(b); j++, i++)
            data[j] = in[i];
    }
}


void NativeThyraVector::copyOutRawData ( double *out ) const
{
    size_t i=0;
    for (size_t b=0; b<numberOfDataBlocks(); b++) {
        const double *data = reinterpret_cast<const double*>(getRawDataBlockAsVoid(b));
        for (size_t j=0; j<sizeOfDataBlock(b); j++, i++)
            out[i] = data[j];
    }
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
    ThyraVectorWrapper* wrapperVector = dynamic_cast<ThyraVectorWrapper*>(ptr);
    if ( wrapperVector!=NULL ) {
        AMP_INSIST(wrapperVector->numVecs()==1,"Not ready for dealing with multiple copies of the vector yet");
        return wrapperVector->getVec(0)->getRawDataBlock<double>(i);
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
    const ThyraVectorWrapper* wrapperVector = dynamic_cast<const ThyraVectorWrapper*>(ptr);
    if ( wrapperVector!=NULL ) {
        AMP_INSIST(wrapperVector->numVecs()==1,"Not ready for dealing with multiple copies of the vector yet");
        return wrapperVector->getVec(0)->getRawDataBlock<double>(i);
    }
    return NULL;
}


size_t NativeThyraVector::numberOfDataBlocks () const 
{ 
    const Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    if ( dynamic_cast<const Thyra::DefaultSpmdVector<double>*>(ptr)!=NULL )
        return 1;
    const ThyraVectorWrapper* wrapperVector = dynamic_cast<const ThyraVectorWrapper*>(ptr);
    if ( wrapperVector!=NULL )
        return wrapperVector->getVec(0)->numberOfDataBlocks();
    AMP_ERROR( "not finished" );
    return 1;
}


size_t NativeThyraVector::sizeOfDataBlock ( size_t i ) const
{
    const Thyra::VectorBase<double>* ptr = d_thyraVec.get();
    const Thyra::DefaultSpmdVector<double>* spmdVector = 
        dynamic_cast<const Thyra::DefaultSpmdVector<double>*>(ptr);
    if ( spmdVector!=NULL )
        return d_local;
    const ThyraVectorWrapper* wrapperVector = dynamic_cast<const ThyraVectorWrapper*>(ptr);
    if ( wrapperVector!=NULL ) {
        AMP_INSIST(wrapperVector->numVecs()==1,"Not ready for dealing with multiple copies of the vector yet");
        return wrapperVector->getVec(0)->sizeOfDataBlock(i);
    }
    AMP_ERROR( "not finished" );
    return d_local;
}


}
}

