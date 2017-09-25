#include "discretization/DOF_ManagerParameters.h"

namespace AMP {
namespace Discretization {


/************************************************************************
 *  Constructors                                                         *
 ************************************************************************/
DOFManagerParameters::DOFManagerParameters() {}
DOFManagerParameters::DOFManagerParameters( AMP::shared_ptr<AMP::Mesh::Mesh> mesh_in )
{
    mesh = mesh_in;
}


} // namespace Discretization
} // namespace AMP
