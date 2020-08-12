#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/trilinos/epetra/ManagedEpetraVectorOperations.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "EpetraVectorEngine.h"


namespace AMP {
namespace LinearAlgebra {


ManagedEpetraVector::ManagedEpetraVector( VectorParameters::shared_ptr params )
    : ManagedVector( params ), EpetraVector()
{
  d_VectorOps = new ManagedEpetraVectorOperations();
}


ManagedEpetraVector::ManagedEpetraVector( shared_ptr alias )
    : ManagedVector( alias ), EpetraVector()
{
  d_VectorOps = new ManagedEpetraVectorOperations();
}

ManagedEpetraVector::~ManagedEpetraVector()
{
  delete d_VectorOps;
}
  
inline ManagedVector *ManagedEpetraVector::getNewRawPtr() const
{
    return new ManagedEpetraVector( std::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}


inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto p   = std::make_shared<ManagedVectorParameters>();
    auto vec = std::dynamic_pointer_cast<Vector>( d_Engine );
    if ( vec ) {
        auto vec2   = vec->cloneVector( "ManagedPetscVectorClone" );
        p->d_Buffer = std::dynamic_pointer_cast<VectorData>( vec2 );
        p->d_Engine = std::dynamic_pointer_cast<Vector>( vec2 );
    } else {
        AMP_ERROR( "ManagedPetscVector::rawClone() should not have reached here!" );
    }
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
