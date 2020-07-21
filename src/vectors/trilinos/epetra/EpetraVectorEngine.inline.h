

namespace AMP {
namespace LinearAlgebra {


inline EpetraVectorEngine::~EpetraVectorEngine() {}

inline Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() { return d_epetraVector; }

inline const Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() const { return d_epetraVector; }

inline AMP_MPI EpetraVectorEngine::getComm() const {
  #if 1
  return VectorOperations::getComm();
#else
  return getEngineParameters()->getComm();
#endif
}

inline bool EpetraVectorEngine::sameEngine( VectorEngine &e ) const
{
    return dynamic_cast<EpetraVectorEngine *>( &e ) != nullptr;
}


} // namespace LinearAlgebra
} // namespace AMP
