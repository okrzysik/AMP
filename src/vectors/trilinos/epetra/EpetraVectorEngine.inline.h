#include "utils/Counter.h"


namespace AMP {
namespace LinearAlgebra {


inline EpetraVectorEngine::~EpetraVectorEngine() { ADD_COUNT( "FLOPS", d_epetraVector.Flops() ); }

inline Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() { return d_epetraVector; }

inline const Epetra_Vector &EpetraVectorEngine::getEpetra_Vector() const { return d_epetraVector; }

inline AMP_MPI EpetraVectorEngine::getComm() const { return getEngineParameters()->getComm(); }

inline void *EpetraVectorEngine::getDataBlock( size_t i )
{
    if ( i > 1 )
        return nullptr;
    double *p;
    getEpetra_Vector().ExtractView( &p );
    return p;
}

inline const void *EpetraVectorEngine::getDataBlock( size_t i ) const
{
    if ( i > 1 )
        return nullptr;
    double *p;
    getEpetra_Vector().ExtractView( &p );
    return p;
}

inline size_t EpetraVectorEngine::getLocalSize() const
{
    return d_Params->castTo<EpetraVectorEngineParameters>().getLocalSize();
}

inline size_t EpetraVectorEngine::getGlobalSize() const
{
    return d_Params->castTo<EpetraVectorEngineParameters>().getGlobalSize();
}

inline bool EpetraVectorEngine::sameEngine( VectorEngine &e ) const
{
    return e.isA<EpetraVectorEngine>();
}

inline size_t EpetraVectorEngine::numberOfDataBlocks() const { return 1; }

inline size_t EpetraVectorEngine::sizeOfDataBlock( size_t i ) const
{
    if ( i != 0 )
        return 0;
    return getLocalSize();
}
}
}
