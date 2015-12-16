#include "vectors/trilinos/EpetraVector.h"
#include "vectors/MultiVector.h"
#include "vectors/trilinos/ManagedEpetraVector.h"


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
    if ( inVector->isA<EpetraVector>() ) {
        retVal = inVector;
    } else if ( inVector->isA<MultiVector>() ) {
        if ( inVector->numberOfDataBlocks() == 1 ) {
            Vector::shared_ptr localVector = inVector->castTo<MultiVector>().getVector( 0 );
            retVal                         = view( localVector );
        } else {
            AMP_ERROR( "View of multi-block MultiVector is not supported yet" );
        }
    } else if ( inVector->isA<ManagedVector>() ) {
        AMP::shared_ptr<Vector> root = inVector->castTo<ManagedVector>().getRootVector();
        if ( root == inVector ) {
            AMP::shared_ptr<ManagedEpetraVector> managed( new ManagedEpetraVector( root ) );
            retVal = managed;
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
    if ( inVector->isA<EpetraVector>() ) {
        return inVector;
    } else if ( inVector->isA<MultiVector>() ) {
        if ( inVector->numberOfDataBlocks() == 1 ) {
            AMP::shared_ptr<MultiVector> multivector = AMP::dynamic_pointer_cast<MultiVector>(
                AMP::const_pointer_cast<Vector>( inVector ) );
            retVal = constView( multivector->getVector( 0 ) );
        } else {
            AMP_ERROR( "View of multi-block MultiVector is not supported yet" );
        }
    } else if ( inVector->isA<ManagedVector>() ) {
        AMP::shared_ptr<ManagedVector> managedVector =
            AMP::dynamic_pointer_cast<ManagedVector>( AMP::const_pointer_cast<Vector>( inVector ) );
        AMP::shared_ptr<Vector> root = managedVector->getRootVector();
        if ( root == inVector ) {
            AMP::shared_ptr<ManagedEpetraVector> managed( new ManagedEpetraVector( root ) );
            retVal = managed;
        } else {
            retVal = constView( root );
        }
    }
    if ( !retVal )
        AMP_ERROR( "Cannot create view!" );
    return retVal;
}
}
}
