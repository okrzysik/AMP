#include "AMP/discretization/DOF_ManagerParameters.h"

namespace AMP {
namespace Discretization {


/************************************************************************
 *  Constructors                                                         *
 ************************************************************************/
DOFManagerParameters::DOFManagerParameters() = default;
DOFManagerParameters::DOFManagerParameters( std::shared_ptr<AMP::Mesh::Mesh> mesh_in )
    : mesh( mesh_in )
{
}


} // namespace Discretization
} // namespace AMP
