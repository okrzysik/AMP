
#include "DTKAMPMeshEntity.h"
#include "DTKAMPMeshEntityImpl.h"

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK Entity interface.
*/

// Constructor.
AMPMeshEntity::AMPMeshEntity( const AMP::Mesh::MeshElement& element,
			      const unsigned long int element_id  )
{
this->b_entity_impl = 
    Teuchos::rcp( new AMPMeshEntityImpl(element, element_id) );
}


}
}
