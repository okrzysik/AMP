#include "ManagedEpetraVector.h"
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
} // namespace LinearAlgebra
} // namespace AMP
