
#ifndef included_AMP_DTK_AMPMeshEntityExtraData
#define included_AMP_DTK_AMPMeshEntityExtraData

#include "AMP/ampmesh/MeshElement.h"

#include "AMP/utils/AMP_MPI.h"

#include <DTK_EntityExtraData.hpp>

namespace AMP {
namespace Operator {


/**
 * AMP Mesh element implementation for DTK EntityExtraData interface.
 */
class AMPMeshEntityExtraData : public DataTransferKit::EntityExtraData
{
public:
    /**
     * Constructor.
     */
    explicit AMPMeshEntityExtraData( const AMP::Mesh::MeshElement &element ) : d_element( element )
    { /* ... */
    }

    //! Destructor
    ~AMPMeshEntityExtraData() {}

    // Underlying mesh element.
    AMP::Mesh::MeshElement d_element;
};
} // namespace Operator
} // namespace AMP

#endif
