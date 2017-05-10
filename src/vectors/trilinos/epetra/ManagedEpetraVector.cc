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

void ManagedEpetraVector::copyVector( Vector::const_shared_ptr vec )
{
    // there must be a more sensible way of doing this but I can't find the documentation - BP
    auto epetraVec = AMP::dynamic_pointer_cast<const ManagedEpetraVector>( vec );
    if ( epetraVec ) {
        double scale = 1.0;
        getEpetra_Vector().Scale( scale, epetraVec->getEpetra_Vector() );
        copyGhostValues( vec );
    } else {
        Vector::copyVector( vec );
    }
}
}
}
