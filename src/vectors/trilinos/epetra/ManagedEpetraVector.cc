#include "vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "vectors/data/VectorDataCPU.h"
#include "utils/Utilities.h"


namespace AMP {
namespace LinearAlgebra {


ManagedEpetraVector::ManagedEpetraVector( VectorParameters::shared_ptr params )
    : ManagedVector( params ), EpetraVector()
{
}


ManagedEpetraVector::ManagedEpetraVector( shared_ptr alias )
    : ManagedVector( alias ), EpetraVector()
{
}

void ManagedEpetraVector::copy( const VectorOperations &src )
{
    // there must be a more sensible way of doing this but I can't find the documentation - BP
    auto epetraVec = dynamic_cast<const ManagedEpetraVector *>( &src );
    if ( epetraVec ) {
        double scale = 1.0;
        getEpetra_Vector().Scale( scale, epetraVec->getEpetra_Vector() );
        copyGhostValues( *dynamic_cast<const VectorData *>( &src ) );
    } else {
        VectorOperationsDefault<double>::copy( src );
    }
}

inline ManagedVector *ManagedEpetraVector::getNewRawPtr() const
{
    return new ManagedEpetraVector( AMP::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}


inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto p = AMP::make_shared<ManagedVectorParameters>();
    p->d_Buffer = AMP::make_shared<VectorDataCPU<double>>(
        d_vBuffer->getLocalStartID(), d_vBuffer->getLocalSize(), d_vBuffer->getGlobalSize() );
    p->d_Engine      = d_pParameters->d_Engine->cloneEngine( p->d_Buffer );
    p->d_CommList    = getCommunicationList();
    p->d_DOFManager  = getDOFManager();
    p->d_CloneEngine = false;
    auto retVal = AMP::make_shared<ManagedEpetraVector>( p );
    retVal->setVariable( var );
    return retVal;
}

inline Epetra_Vector &ManagedEpetraVector::getEpetra_Vector()
{
    auto engine = AMP::dynamic_pointer_cast<EpetraVectorEngine>( d_Engine );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}

inline const Epetra_Vector &ManagedEpetraVector::getEpetra_Vector() const
{
    auto engine = AMP::dynamic_pointer_cast<const EpetraVectorEngine>( d_Engine );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}

inline void ManagedEpetraVector::assemble() {}


} // namespace LinearAlgebra
} // namespace AMP
