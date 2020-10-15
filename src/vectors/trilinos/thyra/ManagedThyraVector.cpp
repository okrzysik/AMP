#include "AMP/vectors/trilinos/thyra/ManagedThyraVector.h"
#include "AMP/vectors/MultiVector.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
ManagedThyraVector::ManagedThyraVector( Vector::shared_ptr vec ) : ManagedVector( vec )
{
    d_thyraVec = Teuchos::RCP<Thyra::VectorBase<double>>(
        new ThyraVectorWrapper( std::vector<Vector::shared_ptr>( 1, vec ) ) );
}
ManagedVector *ManagedThyraVector::getNewRawPtr() const
{
    return new ManagedThyraVector( const_cast<ManagedThyraVector *>( this )->getVectorEngine() );
}


/****************************************************************
 * Destructor                                                    *
 ****************************************************************/
ManagedThyraVector::~ManagedThyraVector() = default;


/****************************************************************
 * Return the vector type                                        *
 ****************************************************************/
std::string ManagedThyraVector::ManagedThyraVector::type() const
{
    std::string retVal = "Managed Thyra Vector";
    retVal += ManagedVector::type();
    return retVal;
}


/****************************************************************
 * Clone the vector                                              *
 ****************************************************************/
Vector::shared_ptr ManagedThyraVector::cloneVector( const Variable::shared_ptr var ) const
{
    auto vec                  = getVectorEngine();
    auto vec2                 = vec->cloneVector( "ManagedThyraVectorClone" );
    auto engine               = std::dynamic_pointer_cast<Vector>( vec2 );
    Vector::shared_ptr retVal = Vector::shared_ptr( new ManagedThyraVector( engine ) );
    retVal->setVariable( var );
    return retVal;
}


/****************************************************************
 * Copy the vector                                               *
 ****************************************************************/
void ManagedThyraVector::copyVector( Vector::const_shared_ptr vec )
{
    Vector::shared_ptr engineVec = std::dynamic_pointer_cast<Vector>( getVectorEngine() );
    engineVec->copyVector( vec );
}

} // namespace LinearAlgebra
} // namespace AMP
