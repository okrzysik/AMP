#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "EpetraVectorEngine.h"


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
    return new ManagedEpetraVector( std::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}


inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto p = std::make_shared<ManagedVectorParameters>();

    // at present all the code I have seen has an implicit
    // assumption that d_vBuffer is non null. In future this should change
    // if we allow for the case that this is null
    p->d_Buffer     = d_vBuffer->cloneData();
    p->d_Engine     = cloneVectorEngine( p->d_Buffer );
    p->d_CommList   = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    auto retVal     = std::make_shared<ManagedEpetraVector>( p );
    retVal->setVariable( var );
    return retVal;
}

inline Epetra_Vector &ManagedEpetraVector::getEpetra_Vector()
{
    auto engine = std::dynamic_pointer_cast<EpetraVectorEngine>( d_Engine );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}

inline const Epetra_Vector &ManagedEpetraVector::getEpetra_Vector() const
{
    auto engine = std::dynamic_pointer_cast<const EpetraVectorEngine>( d_Engine );
    AMP_ASSERT( engine != nullptr );
    return engine->getEpetra_Vector();
}

inline void ManagedEpetraVector::assemble() {}


} // namespace LinearAlgebra
} // namespace AMP
