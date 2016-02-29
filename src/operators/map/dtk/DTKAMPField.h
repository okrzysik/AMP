
#ifndef included_AMP_DTK_AMPField
#define included_AMP_DTK_AMPField

#include <vectors/Vector.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <DTK_Types.hpp>
#include <DTK_Field.hpp>

namespace AMP {
namespace Operator {


/**
 * AMP Vector helpers for DTK.
 */
class DTKAMPField : public DataTransferKit::Field
{
  public:
    /**
     * Constructor.
     */
    explicit DTKAMPField( const AMP::LinearAlgebra::Vector::shared_ptr& amp_vector );

    /*!
     * \brief Get the dimension of the field.
     */
    int dimension() const override;

    /*!
     * \brief Get the locally-owned entity support location ids of the field.
     */
    Teuchos::ArrayView<const DataTransferKit::SupportId> 
	getLocalSupportIds() const override;

    /*!
     * \brief Given a local support id and a dimension, read data from the
     * application field.
     */
    double readFieldData( const DataTransferKit::SupportId support_id,
			  const int dimension ) const override;

    /*!
     * \brief Given a local support id, dimension, and field value, write data
     * into the application field.
     */
    void writeFieldData( const DataTransferKit::SupportId support_id,
			 const int dimension,
			 const double data ) override;

    /*!
     * \brief Finalize a field after writing into it.
     */
    void finalizeAfterWrite() override;

  private:

    // The vector over which the field is defined.
    AMP::LinearAlgebra::Vector::shared_ptr d_amp_vector;

    // The support ids over which the field is constructed.
    Teuchos::Array<DataTransferKit::SupportId> d_support_ids;
};
}
}

#endif
