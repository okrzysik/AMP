

namespace AMP {
namespace LinearAlgebra {


inline EpetraVectorEngine::~EpetraVectorEngine() {}

inline Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() { return d_epetraVector; }

inline const Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() const { return d_epetraVector; }

inline AMP_MPI EpetraVectorEngine::getComm() const {
  return d_Params->getComm();
}

} // namespace LinearAlgebra
} // namespace AMP
