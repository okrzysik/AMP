#include "AMP/vectors/trilinos/epetra/ManagedEpetraVectorOperations.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace AMP {
namespace LinearAlgebra {

//**********************************************************************
// Functions that operate on VectorData objects

void ManagedEpetraVectorOperations::copy( const VectorData &src, VectorData &dst )
{
    // there must be a more sensible way of doing this but I can't find the documentation - BP
    auto srcVec = dynamic_cast<const ManagedEpetraVector *>( &src );
    auto dstVec = dynamic_cast<ManagedEpetraVector *>( &dst );
    if ( srcVec && dstVec ) {
        double scale = 1.0;
        dstVec->getEpetra_Vector().Scale( scale, srcVec->getEpetra_Vector() );
        dst.copyGhostValues( src );
    } else {
        VectorOperationsDefault<double>::copy( src, dst );
    }
    dstVec->dataChanged();
}

} // namespace LinearAlgebra
} // namespace AMP
