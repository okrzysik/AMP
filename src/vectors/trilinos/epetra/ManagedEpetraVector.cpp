#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"


namespace AMP {
namespace LinearAlgebra {


ManagedEpetraVector::ManagedEpetraVector( shared_ptr alias )
    : ManagedVector( alias ), EpetraVector()
{
}

ManagedEpetraVector::~ManagedEpetraVector() {}

inline ManagedVector *ManagedEpetraVector::getNewRawPtr() const
{
    return new ManagedEpetraVector( const_cast<ManagedEpetraVector *>( this )->getVectorEngine() );
}


inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto vec    = getVectorEngine();
    auto vec2   = vec->cloneVector( "ManagedEeptraVectorClone" );
    auto retVal = std::make_shared<ManagedEpetraVector>( vec2 );
    retVal->setVariable( var );
    return retVal;
}

void ManagedEpetraVector::copyVector( Vector::const_shared_ptr vec )
{
    auto engineVec = getVectorEngine();
    engineVec->copyVector( vec );
}

inline Epetra_Vector &ManagedEpetraVector::getEpetra_Vector()
{
    auto data = std::dynamic_pointer_cast<EpetraVectorData>( getVectorEngine()->getVectorData() );
    AMP_ASSERT( data != nullptr );
    return data->getEpetra_Vector();
}

inline const Epetra_Vector &ManagedEpetraVector::getEpetra_Vector() const
{
    auto data =
        std::dynamic_pointer_cast<const EpetraVectorData>( getVectorEngine()->getVectorData() );
    AMP_ASSERT( data != nullptr );
    return data->getEpetra_Vector();
}


} // namespace LinearAlgebra
} // namespace AMP
