#ifndef included_AMP_ThyraVectorSpaceWrapper
#define included_AMP_ThyraVectorSpaceWrapper

// AMP includes
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/trilinos/thyra/ThyraVectorWrapper.h"


// Trilinos includes
DISABLE_WARNINGS
#include "Thyra_VectorSpaceBase.hpp"
ENABLE_WARNINGS


namespace AMP {
namespace LinearAlgebra {


/** \class ThyraVectorSpaceWrapper
 * \brief  Wrapper a VectorSpace in AMP
 * \details  This function is used to allow us to safely wrap an AMP vector
 *   in a thyra vector for use within Trilinos.
 */
class ThyraVectorSpaceWrapper : public Thyra::VectorSpaceBase<double>
{
public:
    //! Default constructor
    explicit ThyraVectorSpaceWrapper( std::shared_ptr<const ThyraVectorWrapper> thyra_vec,
                                      bool is_range = true );

    // Deleted constructors
    ThyraVectorSpaceWrapper()                                  = delete;
    ThyraVectorSpaceWrapper( const ThyraVectorSpaceWrapper & ) = delete;
    ThyraVectorSpaceWrapper &operator=( const ThyraVectorSpaceWrapper & ) = delete;

    //! Destructor
    virtual ~ThyraVectorSpaceWrapper();

    // Virtual functions inherited from VectorSpaceBase
    virtual Teuchos::Ordinal dim() const override;
    virtual bool isCompatible( const Thyra::VectorSpaceBase<double> &vecSpc ) const override;
    virtual Teuchos::RCP<const Thyra::VectorSpaceFactoryBase<double>>
    smallVecSpcFcty() const override;
    virtual double scalarProd( const Thyra::VectorBase<double> &x,
                               const Thyra::VectorBase<double> &y ) const override;


protected:
    // Virtual functions inherited from VectorSpaceBase
    virtual Teuchos::RCP<Thyra::VectorBase<double>> createMember() const override;
    virtual Teuchos::RCP<Thyra::MultiVectorBase<double>>
    createMembers( int numMembers ) const override;
    virtual Teuchos::RCP<Thyra::VectorBase<double>>
    createMemberView( const RTOpPack::SubVectorView<double> &raw_v ) const override;
    virtual Teuchos::RCP<const Thyra::VectorBase<double>>
    createMemberView( const RTOpPack::ConstSubVectorView<double> &raw_v ) const override;
    virtual Teuchos::RCP<Thyra::MultiVectorBase<double>>
    createMembersView( const RTOpPack::SubMultiVectorView<double> &raw_mv ) const override;
    virtual Teuchos::RCP<const Thyra::MultiVectorBase<double>>
    createMembersView( const RTOpPack::ConstSubMultiVectorView<double> &raw_mv ) const override;
    virtual void scalarProdsImpl( const Thyra::MultiVectorBase<double> &X,
                                  const Thyra::MultiVectorBase<double> &Y,
                                  const Teuchos::ArrayView<double> &scalarProds ) const override;

    // Local data
    bool d_is_range;
    std::shared_ptr<const ThyraVectorWrapper> d_thyra_vec;
};


} // namespace LinearAlgebra
} // namespace AMP

#endif
