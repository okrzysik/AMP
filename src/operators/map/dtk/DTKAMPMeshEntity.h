
#ifndef included_AMP_DTK_AMPMeshEntity
#define included_AMP_DTK_AMPMeshEntity

#include "utils/shared_ptr.h"
#include "matrices/Matrix.h"
#include "operators/Operator.h"
#include "operators/OperatorParameters.h"
#include "vectors/Vector.h"

#include <DTK_Entity.hpp>

namespace AMP {
namespace Operator {


/**
  * AMP Mesh element implementation for DTK EntityExtraData interface.
*/
class AMPMeshEntity : public DataTransferKit::Entity
{
public :

    /**
     * Constructor.
     */
    AMPMeshEntity( const AMP::Mesh::MeshElement& element,
		   const unsigned long int element_id );

    //! Destructor
    ~AMPMeshEntity() { }
};


}
}

#endif


