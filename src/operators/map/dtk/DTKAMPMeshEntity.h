
#ifndef included_AMP_DTK_AMPMeshEntity
#define included_AMP_DTK_AMPMeshEntity

#include "ampmesh/MeshElement.h"

#include "utils/AMP_MPI.h"

#include <DTK_Entity.hpp>

#include <unordered_map>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK Entity interface.
*/
class AMPMeshEntity : public DataTransferKit::Entity
{
public:
    /**
     * Constructor.
     */
    explicit AMPMeshEntity( const AMP::Mesh::MeshElement &element,
			    const std::unordered_map<int,int>& rank_map );

    //! Destructor
    ~AMPMeshEntity() {}
};
}
}

#endif
