#include "AMP/vectors/trilinos/thyra/ThyraVector.h"

#include "AMP/vectors/MultiVariable.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/trilinos/thyra/ManagedThyraVector.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
ThyraVector::ThyraVector() { d_thyraVec.reset(); }


/****************************************************************
 * Destructors                                                   *
 ****************************************************************/
ThyraVector::~ThyraVector() { d_thyraVec.reset(); }


/****************************************************************
 * view                                                          *
 ****************************************************************/
Vector::const_shared_ptr ThyraVector::constView( Vector::const_shared_ptr inVector )
{
    return view( std::const_pointer_cast<Vector>( inVector ) );
}
Vector::shared_ptr ThyraVector::view( Vector::shared_ptr inVector )
{
    // Check if we have an exisiting view
    if ( std::dynamic_pointer_cast<ThyraVector>( inVector ) != nullptr )
        return inVector;
    if ( inVector->hasView<ManagedThyraVector>() )
        return inVector->getView<ManagedThyraVector>();
    // Create a new view
    Vector::shared_ptr retVal;
    if ( std::dynamic_pointer_cast<ManagedVector>( inVector ) ) {
        retVal = std::make_shared<ManagedThyraVector>( inVector );
        inVector->registerView( retVal );
    } else if ( std::dynamic_pointer_cast<MultiVector>( inVector ) ) {
        auto newParams      = std::make_shared<ManagedThyraVectorParameters>();
        newParams->d_Engine = std::dynamic_pointer_cast<Vector>( inVector );
        AMP_INSIST( inVector->getCommunicationList().get() != nullptr,
                    "All vectors must have a communication list" );
        newParams->d_CommList = inVector->getCommunicationList();
        AMP_INSIST( inVector->getDOFManager().get() != nullptr,
                    "All vectors must have a DOFManager list" );
        newParams->d_DOFManager = inVector->getDOFManager();
        auto newVector          = std::make_shared<ManagedThyraVector>( newParams );
        newVector->setVariable( inVector->getVariable() );
        newVector->getVectorData()->setUpdateStatusPtr(
            inVector->getVectorData()->getUpdateStatusPtr() );
        retVal = newVector;
        inVector->registerView( retVal );
    } else {
        retVal = view( MultiVector::view( inVector, inVector->getComm() ) );
        inVector->registerView( retVal );
    }
    return retVal;
}


/****************************************************************
 * Return the thyra vector                                       *
 ****************************************************************/
Teuchos::RCP<Thyra::VectorBase<double>> ThyraVector::getVec() { return d_thyraVec; }
Teuchos::RCP<const Thyra::VectorBase<double>> ThyraVector::getVec() const { return d_thyraVec; }


/****************************************************************
 * Return the views to the AMP vectors                           *
 ****************************************************************/
template<class T>
static void nullDeleter( T * )
{
}
AMP::LinearAlgebra::Vector::shared_ptr ThyraVector::view( Thyra::VectorBase<double> *vec )
{
    AMP::LinearAlgebra::Vector::shared_ptr vec_out;
    if ( vec == nullptr ) {
        // Null vec, do nothing
    } else if ( dynamic_cast<AMP::LinearAlgebra::ThyraVectorWrapper *>( vec ) ) {
        auto *tmp = dynamic_cast<AMP::LinearAlgebra::ThyraVectorWrapper *>( vec );
        if ( tmp->numVecs() == 0 ) {
            vec_out.reset();
        } else if ( tmp->numVecs() == 1 ) {
            vec_out = tmp->getVec( 0 );
        } else {
            std::vector<AMP::LinearAlgebra::Variable::shared_ptr> vars;
            for ( size_t i = 0; i < tmp->d_vecs.size(); i++ ) {
                char name[100];
                sprintf( name, "col-%i\n", (int) tmp->d_cols[i] );
                vars.push_back( std::make_shared<AMP::LinearAlgebra::Variable>( name ) );
            }
            auto multiVar =
                std::make_shared<AMP::LinearAlgebra::MultiVariable>( "ThyraMultiVec", vars );
            vec_out = AMP::LinearAlgebra::MultiVector::create(
                multiVar, tmp->d_vecs[0]->getComm(), tmp->d_vecs );
            // Currently our multivectors can't be easily subsetted to create the original vectors
            AMP_ERROR( "Not ready for ThyraMultiVectors yet" );
        }
    } else {
        AMP_ERROR( "Not finished" );
    }
    return vec_out;
}
AMP::LinearAlgebra::Vector::const_shared_ptr
ThyraVector::constView( const Thyra::VectorBase<double> *vec )
{
    AMP::LinearAlgebra::Vector::const_shared_ptr vec_out;
    if ( vec == nullptr ) {
        // Null vec, do nothing
    } else if ( dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper *>( vec ) ) {
        const auto *tmp = dynamic_cast<const AMP::LinearAlgebra::ThyraVectorWrapper *>( vec );
        if ( tmp->numVecs() == 0 ) {
            vec_out.reset();
        } else if ( tmp->numVecs() == 1 ) {
            vec_out = tmp->getVec( 0 );
        } else {
            std::vector<AMP::LinearAlgebra::Variable::shared_ptr> vars;
            for ( size_t i = 0; i < tmp->d_vecs.size(); i++ ) {
                char name[100];
                sprintf( name, "col-%i\n", (int) tmp->d_cols[i] );
                vars.push_back( std::make_shared<AMP::LinearAlgebra::Variable>( name ) );
            }
            auto multiVar =
                std::make_shared<AMP::LinearAlgebra::MultiVariable>( "ThyraMultiVec", vars );
            vec_out = AMP::LinearAlgebra::MultiVector::create(
                multiVar, tmp->d_vecs[0]->getComm(), tmp->d_vecs );
            // Currently our multivectors can't be easily subsetted to create the original vectors
            AMP_ERROR( "Not ready for ThyraMultiVectors yet" );
        }
    } else {
        AMP_ERROR( "Not finished" );
    }
    return vec_out;
}


} // namespace LinearAlgebra
} // namespace AMP
