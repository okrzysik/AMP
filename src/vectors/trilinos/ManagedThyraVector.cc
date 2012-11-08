#include "vectors/trilinos/ManagedThyraVector.h"
#include "vectors/SimpleVector.h"


namespace AMP {
namespace LinearAlgebra {



/****************************************************************
* Constructors                                                  *
****************************************************************/
ManagedThyraVector::ManagedThyraVector( VectorParameters::shared_ptr  params ):
    ManagedVector( params )
{
    AMP_ERROR("Not implimented yet");
}
ManagedThyraVector::ManagedThyraVector( Vector::shared_ptr  alias ):
    ManagedVector( alias )
{
    AMP_ERROR("Not implimented yet");
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
    /*boost::shared_ptr<ManagedVectorParameters>  p ( new ManagedVectorParameters () );
    p->d_Buffer = VectorEngine::BufferPtr ( new VectorEngine::Buffer ( d_vBuffer->size() ) );
    p->d_Engine = d_pParameters->d_Engine->cloneEngine( p->d_Buffer );
    p->d_CommList = getCommunicationList();
    p->d_DOFManager = getDOFManager();
    p->d_CloneEngine = false;
    Vector::shared_ptr retVal = Vector::shared_ptr ( new ManagedEpetraVector ( boost::dynamic_pointer_cast<VectorParameters> ( p ) ) );
    retVal->setVariable ( var );
    return retVal;*/
    AMP_ERROR("Not implimented yet");
    return Vector::shared_ptr();
}


/****************************************************************
* Copy the vector                                               *
****************************************************************/
void ManagedThyraVector::copyVector(const Vector::const_shared_ptr &vec)
{
    AMP_ERROR("Not implimented yet");
}



}
}

