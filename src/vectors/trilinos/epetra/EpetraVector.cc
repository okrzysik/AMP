#include "vectors/trilinos/epetra/EpetraVector.h"
#include "vectors/MultiVector.h"
#include "vectors/trilinos/epetra/ManagedEpetraVector.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
* Constructors / De-constructors                        *
********************************************************/
EpetraVector::EpetraVector() {}
EpetraVector::~EpetraVector() {}


/********************************************************
* View                                                  *
********************************************************/
Vector::shared_ptr EpetraVector::view( Vector::shared_ptr inVector )
{
    Vector::shared_ptr retVal;
    if ( dynamic_pointer_cast<EpetraVector>(inVector) ) {
        retVal = inVector;
    } else if ( dynamic_pointer_cast<MultiVector>(inVector) ) {
        auto multivec = dynamic_pointer_cast<MultiVector>(inVector);
        if ( inVector->numberOfDataBlocks() == 1 ) {
            retVal = view( multivec->getVector( 0 ) );
        } else {
            AMP_ERROR( "View of multi-block MultiVector is not supported yet" );
        }
    } else if ( dynamic_pointer_cast<ManagedVector>(inVector) ) {
        auto managed = dynamic_pointer_cast<ManagedVector>(inVector);
        auto root = managed->getRootVector();
        if ( root == inVector ) {
            retVal = AMP::make_shared<ManagedEpetraVector>( root );
        } else {
            retVal = view( root );
        }
    }
    if ( !retVal )
        AMP_ERROR( "Cannot create view!" );
    return retVal;
}
Vector::const_shared_ptr EpetraVector::constView( Vector::const_shared_ptr inVector )
{
    Vector::const_shared_ptr retVal;
    if ( dynamic_pointer_cast<const EpetraVector>(inVector) ) {
        return inVector;
    } else if ( dynamic_pointer_cast<const MultiVector>(inVector) ) {
        auto multivec = dynamic_pointer_cast<const MultiVector>(inVector);
        if ( inVector->numberOfDataBlocks() == 1 ) {
            retVal = constView( multivec->getVector( 0 ) );
        } else {
            AMP_ERROR( "View of multi-block MultiVector is not supported yet" );
        }
    } else if ( dynamic_pointer_cast<const ManagedVector>(inVector) ) {
        auto managed = dynamic_pointer_cast<const ManagedVector>(inVector);
        auto root = AMP::const_pointer_cast<ManagedVector>(managed)->getRootVector();
        if ( root == inVector ) {
            retVal = AMP::make_shared<ManagedEpetraVector>( root );
        } else {
            retVal = constView( root );
        }
    } else {
        AMP::shared_ptr<ManagedEpetraVector> managed( new ManagedEpetraVector( AMP::const_pointer_cast<Vector>( inVector ) ) );
        retVal = managed;
    }

    if ( !retVal )
        AMP_ERROR( "Cannot create view!" );
    return retVal;
}
}
}
