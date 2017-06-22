namespace AMP {
namespace LinearAlgebra {

inline ManagedVector *ManagedEpetraVector::getNewRawPtr() const
{
    return new ManagedEpetraVector( AMP::dynamic_pointer_cast<VectorParameters>( d_pParameters ) );
}

inline std::string ManagedEpetraVector::type() const
{
    std::string retVal = "Managed Epetra Vector";
    retVal += ManagedVector::type();
    return retVal;
}
inline Vector::shared_ptr ManagedEpetraVector::cloneVector( const Variable::shared_ptr var ) const
{
    AMP::shared_ptr<ManagedVectorParameters> p( new ManagedVectorParameters() );
    p->d_Buffer      = VectorEngine::BufferPtr( new VectorEngine::Buffer( d_vBuffer->size() ) );
    p->d_Engine      = d_pParameters->d_Engine->cloneEngine( p->d_Buffer );
    p->d_CommList    = getCommunicationList();
    p->d_DOFManager  = getDOFManager();
    p->d_CloneEngine = false;
    Vector::shared_ptr retVal = Vector::shared_ptr(
        new ManagedEpetraVector( AMP::dynamic_pointer_cast<VectorParameters>( p ) ) );
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
}
}
