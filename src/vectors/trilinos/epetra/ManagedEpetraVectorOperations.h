#ifndef included_AMP_ManagedEpetraVectorOperations
#define included_AMP_ManagedEpetraVectorOperations

#include "AMP/vectors/operations/ManagedVectorOperations.h"


namespace AMP {
namespace LinearAlgebra {

/**
   \brief Class used to control data and kernels of various vector libraries
   \details  A ManagedVector will take an engine and create a buffer, if
   necessary.

   A ManagedVector has two pointers: data and engine.  If the data pointer
   is null, then the engine is assumed to have the data.
*/
class ManagedEpetraVectorOperations : public ManagedVectorOperations
{

public:
    ManagedEpetraVectorOperations() : ManagedVectorOperations(){};

public:
    //**********************************************************************
    // functions that operate on VectorData
    void copy( const VectorData &src, VectorData &dst ) override;

public: // Pull VectorOperations into the current scope
    using ManagedVectorOperations::copy;
};


} // namespace LinearAlgebra
} // namespace AMP


#endif
