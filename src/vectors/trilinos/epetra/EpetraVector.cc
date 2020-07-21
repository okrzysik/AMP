#include "AMP/vectors/trilinos/epetra/EpetraVector.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/utils/Utilities.h"


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
createManagedEpetraVector( Vector::shared_ptr inVector, std::shared_ptr<EpetraVectorEngine> engine )
{
    auto newParams           = std::make_shared<ManagedVectorParameters>();
    newParams->d_Engine      = engine;
    newParams->d_CloneEngine = false;
    AMP_INSIST( inVector->getCommunicationList(), "All vectors must have a communication list" );
    newParams->d_CommList = inVector->getCommunicationList();
    AMP_INSIST( inVector->getDOFManager(), "All vectors must have a DOFManager list" );
    newParams->d_DOFManager = inVector->getDOFManager();
    auto retVal             = std::make_shared<ManagedEpetraVector>( newParams );
    retVal->setVariable( inVector->getVariable() );
    retVal->setUpdateStatusPtr( inVector->getUpdateStatusPtr() );
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
    } else if ( std::dynamic_pointer_cast<EpetraVectorEngine>( inVector ) ) {
        auto engine = std::dynamic_pointer_cast<EpetraVectorEngine>( inVector );
        retVal      = createManagedEpetraVector( inVector, engine );
    } else {
        // Create a multivector to wrap the given vector and create a view
        auto engineParams = std::make_shared<EpetraVectorEngineParameters>(
            inVector->getLocalSize(), inVector->getGlobalSize(), inVector->getComm() );
        auto engine = std::make_shared<EpetraVectorEngine>( engineParams, inVector );
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
