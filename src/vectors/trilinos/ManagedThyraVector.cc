#include "vectors/trilinos/ManagedThyraVector.h"
#include "vectors/SimpleVector.h"
#include "vectors/MultiVector.h"

#include "vectors/trilinos/ThyraVectorWrapper.h"


namespace AMP {
namespace LinearAlgebra {



/****************************************************************
* Constructors                                                  *
****************************************************************/
ManagedThyraVector::ManagedThyraVector( VectorParameters::shared_ptr  params ):
    ManagedVector( params )
{
    Vector::shared_ptr vec = boost::dynamic_pointer_cast<Vector>( d_Engine );
    d_thyraVec = Teuchos::RCP<Thyra::VectorBase<double> >( new ThyraVectorWrapper(vec) );
}
ManagedThyraVector::ManagedThyraVector( Vector::shared_ptr  alias ):
    ManagedVector( alias )
{
    d_thyraVec = Teuchos::RCP<Thyra::VectorBase<double> >( new ThyraVectorWrapper(alias) );
}
ManagedVector* ManagedThyraVector::getNewRawPtr () const
{ 
    return new ManagedThyraVector( boost::dynamic_pointer_cast<VectorParameters>( d_pParameters ) ); 
}


/****************************************************************
* Destructor                                                    *
****************************************************************/
ManagedThyraVector::~ManagedThyraVector( )
{
}


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
Vector::shared_ptr  ManagedThyraVector::cloneVector ( const Variable::shared_ptr var ) const
{
    boost::shared_ptr<ManagedThyraVectorParameters>  p ( new ManagedThyraVectorParameters() );
    p->d_Engine = d_pParameters->d_Engine->cloneEngine( p->d_Buffer );
    p->d_CommList = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    p->d_CloneEngine = false;
    Vector::shared_ptr retVal = Vector::shared_ptr( new ManagedThyraVector( p ) );
    retVal->setVariable ( var );
    return retVal;
}


/****************************************************************
* Copy the vector                                               *
****************************************************************/
void ManagedThyraVector::copyVector(const Vector::const_shared_ptr &vec)
{
    Vector::shared_ptr engineVec = boost::dynamic_pointer_cast<Vector>( d_Engine );
    engineVec->copyVector( vec );
}



}
}

