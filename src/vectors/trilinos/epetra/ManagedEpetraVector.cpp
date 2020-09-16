#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "EpetraVectorEngine.h"


namespace AMP {
namespace LinearAlgebra {


ManagedEpetraVector::ManagedEpetraVector( std::shared_ptr<ManagedVectorParameters> params )
    : ManagedVector( params ), EpetraVector()
{
}


ManagedEpetraVector::ManagedEpetraVector( shared_ptr alias )
    : ManagedVector( alias ), EpetraVector()
{
}

ManagedEpetraVector::~ManagedEpetraVector() {}

inline ManagedVector *ManagedEpetraVector::getNewRawPtr() const
{
    return new ManagedEpetraVector( d_pParameters );
}


inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto p   = std::make_shared<ManagedVectorParameters>();
    auto vec = getVectorEngine();
    if ( vec ) {
        auto vec2   = vec->cloneVector( "ManagedEeptraVectorClone" );
        p->d_Buffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        p->d_Engine = std::dynamic_pointer_cast<Vector>( vec2 );
    } else {
        AMP_ERROR( "ManagedEpetraVector::rawClone() should not have reached here!" );
    }
    p->d_CommList   = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    auto retVal     = std::make_shared<ManagedEpetraVector>( p );
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
    auto engine = std::dynamic_pointer_cast<EpetraVectorEngine>( getVectorEngine() );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}

inline const Epetra_Vector &ManagedEpetraVector::getEpetra_Vector() const
{
    auto engine = std::dynamic_pointer_cast<const EpetraVectorEngine>( getVectorEngine() );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}


} // namespace LinearAlgebra
} // namespace AMP
