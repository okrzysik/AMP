
#include "DTKAMPField.h"

namespace AMP {
namespace Operator {

//---------------------------------------------------------------------------//
// Constructor.
DTKAMPField::DTKAMPField( const AMP::shared_ptr<AMP::LinearAlgebra::Vector> &amp_vector )
    : d_amp_vector( amp_vector ), d_support_ids( 0 )
{
    // Get the dof IDs if amp_vector is not NULL
    if ( d_amp_vector.get() != NULL ) {
        // Get the dof manager
        AMP::shared_ptr<AMP::Discretization::DOFManager> dof_manager;
        dof_manager = d_amp_vector->getDOFManager();

        // Get the number of degrees of freedom
        AMP::Mesh::MeshIterator meshIterator = dof_manager->getIterator();
        const std::size_t n                  = meshIterator.size();
        d_support_ids.reserve( n );

        // Iterate over the DOF manager elements and extract their DOF ids.
        std::vector<std::size_t> dofIndices;
        meshIterator = meshIterator.begin();
        for ( std::size_t i = 0; i < n; ++i, ++meshIterator ) {
            if ( meshIterator->globalID().is_local() ) {
                dof_manager->getDOFs( meshIterator->globalID(), dofIndices );
                AMP_ASSERT( dofIndices.size() == 1 );
                d_support_ids.push_back( dofIndices[0] );
            }
        }
        AMP_ASSERT( meshIterator == meshIterator.end() );
    }
}

//---------------------------------------------------------------------------//
// Get the dimension of the field.
int DTKAMPField::dimension() const
{
    // Only 1D suppoprted for now.
    return 1;
}

//---------------------------------------------------------------------------//
// Get the locally-owned entity support ids of the field.
Teuchos::ArrayView<const DataTransferKit::SupportId> DTKAMPField::getLocalSupportIds() const
{
    return d_support_ids();
}

//---------------------------------------------------------------------------//
// Given a support id and dimension, read data from the field.
double DTKAMPField::readFieldData( const DataTransferKit::SupportId support_id,
                                   const int dimension ) const
{
    AMP_ASSERT( 0 == dimension );
    return d_amp_vector->getValueByGlobalID( support_id );
}

//---------------------------------------------------------------------------//
// Given an support id, dimension, and field value, write the field vaue.
void DTKAMPField::writeFieldData( const DataTransferKit::SupportId support_id,
                                  const int dimension,
                                  const double data )
{
    AMP_ASSERT( 0 == dimension );
    d_amp_vector->setValueByGlobalID( support_id, data );
}

//---------------------------------------------------------------------------//
// Finalize a field after writing into it.
void DTKAMPField::finalizeAfterWrite()
{
    // Do nothing for now.
    if ( d_amp_vector )
        d_amp_vector->makeConsistent( AMP::LinearAlgebra::Vector::ScatterType::CONSISTENT_SET );
}

//---------------------------------------------------------------------------//

} // end namespace Operator
} // end namespace AMP
