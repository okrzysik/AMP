#include "AMP/vectors/trilinos/epetra/EpetraVector.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/VectorBuilder.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"
#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"


namespace AMP {
namespace LinearAlgebra {


/********************************************************
 * Constructors / De-constructors                        *
 ********************************************************/
EpetraVector::EpetraVector()  = default;
EpetraVector::~EpetraVector() = default;


/********************************************************
 * View                                                  *
 ********************************************************/
static std::shared_ptr<ManagedEpetraVector>
createManagedEpetraVector( Vector::shared_ptr inVector, std::shared_ptr<Vector> engine )
{
    auto retVal = std::make_shared<ManagedEpetraVector>( engine );
    retVal->setVariable( inVector->getVariable() );
    retVal->getVectorData()->setUpdateStatusPtr( inVector->getVectorData()->getUpdateStatusPtr() );
    inVector->registerView( retVal );
    return retVal;
}
Vector::shared_ptr EpetraVector::view( Vector::shared_ptr inVector )
{
    AMP_INSIST( inVector->numberOfDataBlocks() == 1,
                "Epetra does not support more than 1 data block" );
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<EpetraVector>( inVector ) ) {
        retVal = inVector;
    } else if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto multivec = std::dynamic_pointer_cast<MultiVector>( inVector );
        if ( multivec->getNumberOfSubvectors() == 1 ) {
            retVal = view( multivec->getVector( 0 ) );
        } else {
            AMP_ERROR( "View of multi-block MultiVector is not supported yet" );
        }
    } else if ( std::dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        auto managed = std::dynamic_pointer_cast<ManagedVector>( inVector );
        auto root    = managed->getRootVector();
        if ( root == inVector ) {
            retVal = std::make_shared<ManagedEpetraVector>( root );
        } else {
            retVal = view( root );
        }
    } else if ( std::dynamic_pointer_cast<EpetraVectorData>( inVector->getVectorData() ) ) {
        retVal = createManagedEpetraVector( inVector, inVector );
    } else {
        auto engine = createEpetraVector( inVector->getCommunicationList(),
                                          inVector->getDOFManager(),
                                          inVector->getVectorData() );
        retVal      = createManagedEpetraVector( inVector, engine );
    }
    if ( !retVal )
        AMP_ERROR( "Cannot create view!" );
    return retVal;
}
Vector::const_shared_ptr EpetraVector::constView( Vector::const_shared_ptr inVector )
{
    return view( std::const_pointer_cast<Vector>( inVector ) );
}


} // namespace LinearAlgebra
} // namespace AMP
