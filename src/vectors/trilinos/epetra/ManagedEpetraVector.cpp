#include "AMP/vectors/trilinos/epetra/ManagedEpetraVector.h"
#include "AMP/vectors/data/ManagedVectorData.h"
#include "AMP/vectors/data/VectorDataCPU.h"
#include "AMP/vectors/operations/ManagedVectorOperations.h"
#include "AMP/vectors/trilinos/epetra/EpetraVectorData.h"


namespace AMP {
namespace LinearAlgebra {


static inline auto getVectorEngine( const std::shared_ptr<VectorData> &data )
{
    auto managed = std::dynamic_pointer_cast<ManagedVectorData>( data );
    AMP_ASSERT( managed );
    return managed->getVectorEngine();
}
static inline auto getVectorEngine( const std::shared_ptr<const VectorData> &data )
{
    auto managed = std::dynamic_pointer_cast<const ManagedVectorData>( data );
    AMP_ASSERT( managed );
    return managed->getVectorEngine();
}


ManagedEpetraVector::ManagedEpetraVector( shared_ptr vec ) : Vector(), EpetraVector()
{
    AMP_ASSERT( !std::dynamic_pointer_cast<ManagedVectorData>( vec->getVectorData() ) );
    d_VectorOps  = std::make_shared<ManagedVectorOperations>();
    d_VectorData = std::make_shared<ManagedVectorData>( vec );
    d_DOFManager = vec->getDOFManager();
    setVariable( vec->getVariable() );
}

ManagedEpetraVector::~ManagedEpetraVector() {}

std::unique_ptr<Vector> ManagedEpetraVector::rawClone( const Variable::shared_ptr var ) const
{
    auto vec    = getVectorEngine( getVectorData() );
    auto vec2   = vec->cloneVector( "ManagedEeptraVectorClone" );
    auto retVal = std::make_unique<ManagedEpetraVector>( vec2 );
    retVal->setVariable( var );
    return retVal;
}

void ManagedEpetraVector::copyVector( Vector::const_shared_ptr vec )
{
    auto engineVec = getVectorEngine( getVectorData() );
    engineVec->copyVector( vec );
}

Epetra_Vector &ManagedEpetraVector::getEpetra_Vector()
{
    auto vec  = getVectorEngine( getVectorData() );
    auto data = std::dynamic_pointer_cast<EpetraVectorData>( vec->getVectorData() );
    AMP_ASSERT( data != nullptr );
    return data->getEpetra_Vector();
}

const Epetra_Vector &ManagedEpetraVector::getEpetra_Vector() const
{
    auto vec  = getVectorEngine( getVectorData() );
    auto data = std::dynamic_pointer_cast<const EpetraVectorData>( vec->getVectorData() );
    AMP_ASSERT( data != nullptr );
    return data->getEpetra_Vector();
}

void ManagedEpetraVector::swapVectors( Vector &other )
{
    d_VectorData->swapData( *other.getVectorData() );
}


/********************************************************
 * Subset                                                *
 ********************************************************/
Vector::shared_ptr ManagedEpetraVector::subsetVectorForVariable( Variable::const_shared_ptr name )
{
    Vector::shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::subsetVectorForVariable( name );
    if ( !retVal ) {
        auto vec = getVectorEngine( getVectorData() );
        if ( vec )
            retVal = vec->subsetVectorForVariable( name );
    }
    return retVal;
}
Vector::const_shared_ptr
ManagedEpetraVector::constSubsetVectorForVariable( Variable::const_shared_ptr name ) const
{
    Vector::const_shared_ptr retVal;
    if ( !retVal )
        retVal = Vector::constSubsetVectorForVariable( name );
    if ( !retVal ) {
        auto const vec = getVectorEngine( getVectorData() );
        if ( vec )
            retVal = vec->constSubsetVectorForVariable( name );
    }
    if ( !retVal ) {
        auto const vec = getVectorEngine( getVectorData() );
        printf( "Unable to subset for %s in %s:%s\n",
                name->getName().data(),
                getVariable()->getName().data(),
                vec->getVariable()->getName().data() );
    }
    return retVal;
}


} // namespace LinearAlgebra
} // namespace AMP
