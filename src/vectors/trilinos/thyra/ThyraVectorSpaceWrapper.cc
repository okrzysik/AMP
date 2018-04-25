#include "AMP/vectors/trilinos/thyra/ThyraVectorSpaceWrapper.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"

#include "AMP/utils/Utilities.h"

DISABLE_WARNINGS
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_VectorSpaceBase.hpp"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {


/****************************************************************
 * Constructors                                                  *
 ****************************************************************/
ThyraVectorSpaceWrapper::ThyraVectorSpaceWrapper(
    AMP::shared_ptr<const ThyraVectorWrapper> thyra_vec, bool is_range )
{
    AMP_INSIST( thyra_vec != nullptr, "thyra_vec may not be NULL" );
    d_thyra_vec = thyra_vec;
    d_is_range  = is_range;
}


/****************************************************************
 * Destructor                                                    *
 ****************************************************************/
ThyraVectorSpaceWrapper::~ThyraVectorSpaceWrapper() {}


/****************************************************************
 * Virtual functions inherited from VectorSpaceBase              *
 ****************************************************************/
Teuchos::Ordinal ThyraVectorSpaceWrapper::dim() const
{
    if ( !d_is_range )
        return static_cast<Teuchos::Ordinal>( d_thyra_vec->numColumns() );
    return d_thyra_vec->numRows();
}
bool ThyraVectorSpaceWrapper::isCompatible( const Thyra::VectorSpaceBase<double> &vecSpc ) const
{
    const auto *vecSpaceWrapper = dynamic_cast<const ThyraVectorSpaceWrapper *>( &vecSpc );
    if ( vecSpaceWrapper == nullptr )
        return false;
    if ( this == vecSpaceWrapper )
        return true;
    auto dofs1 = d_thyra_vec->getDOFManager();
    auto dofs2 = vecSpaceWrapper->d_thyra_vec->getDOFManager();
    if ( dofs1 == dofs2 )
        return true;
    if ( *dofs1 == *dofs2 )
        return true;
    return false;
}
Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<double>>
ThyraVectorSpaceWrapper::smallVecSpcFcty() const
{
    AMP_ASSERT( d_is_range );
    AMP_ERROR( "Not finished" );
    return Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<double>>();
}
double ThyraVectorSpaceWrapper::scalarProd( const Thyra::VectorBase<double> &x,
                                            const Thyra::VectorBase<double> &y ) const
{
    return dot( x, y );
}
Teuchos::RCP<Thyra::VectorBase<double>> ThyraVectorSpaceWrapper::createMember() const
{
    AMP_ASSERT( d_is_range );
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs( 1 );
    vecs[0] = d_thyra_vec->getVec( 0 )->cloneVector();
    return Teuchos::RCP<Thyra::VectorBase<double>>( new ThyraVectorWrapper( vecs ) );
}
Teuchos::RCP<Thyra::MultiVectorBase<double>>
ThyraVectorSpaceWrapper::createMembers( int numMembers ) const
{
    AMP_ASSERT( d_is_range );
    std::vector<AMP::LinearAlgebra::Vector::shared_ptr> vecs( numMembers );
    for ( int i = 0; i < numMembers; i++ )
        vecs[i] = d_thyra_vec->getVec( 0 )->cloneVector();
    return Teuchos::RCP<Thyra::VectorBase<double>>( new ThyraVectorWrapper( vecs ) );
}
Teuchos::RCP<Thyra::VectorBase<double>>
ThyraVectorSpaceWrapper::createMemberView( const RTOpPack::SubVectorView<double> &raw_v ) const
{
    NULL_USE( raw_v );
    AMP_ERROR( "Not finished" );
    return Teuchos::RCP<Thyra::VectorBase<double>>();
}
Teuchos::RCP<const Thyra::VectorBase<double>>
ThyraVectorSpaceWrapper::createMemberView( const RTOpPack::ConstSubVectorView<double> &raw_v ) const
{
    NULL_USE( raw_v );
    AMP_ERROR( "Not finished" );
    return Teuchos::RCP<const Thyra::VectorBase<double>>();
}
Teuchos::RCP<Thyra::MultiVectorBase<double>> ThyraVectorSpaceWrapper::createMembersView(
    const RTOpPack::SubMultiVectorView<double> &raw_mv ) const
{
    AMP_ASSERT( !d_is_range );
    size_t N_rows = d_thyra_vec->numColumns();
    auto space    = Thyra::defaultSpmdVectorSpace<double>( N_rows );
    auto view     = Thyra::createMembersView<double>( space, raw_mv, "" );
    return view;
}
Teuchos::RCP<const Thyra::MultiVectorBase<double>> ThyraVectorSpaceWrapper::createMembersView(
    const RTOpPack::ConstSubMultiVectorView<double> &raw_mv ) const
{
    AMP_ASSERT( !d_is_range );
    size_t N_rows = d_thyra_vec->numColumns();
    auto space    = Thyra::defaultSpmdVectorSpace<double>( N_rows );
    auto view     = Thyra::createMembersView<double>( space, raw_mv, "" );
    return view;
}
void ThyraVectorSpaceWrapper::scalarProdsImpl( const Thyra::MultiVectorBase<double> &X,
                                               const Thyra::MultiVectorBase<double> &Y,
                                               const Teuchos::ArrayView<double> &scalarProds ) const
{
    NULL_USE( X );
    NULL_USE( Y );
    NULL_USE( scalarProds );
    AMP_ERROR( "Not finished" );
}
} // namespace LinearAlgebra
} // namespace AMP
