#include "discretization/DOF_ManagerParameters.h"

namespace AMP {
namespace Discretization {


/************************************************************************
*  Constructors                                                         *
************************************************************************/
DOFManagerParameters::DOFManagerParameters() 
{
}
DOFManagerParameters::DOFManagerParameters( boost::shared_ptr<AMP::Mesh::Mesh> mesh_in )
{
    mesh = mesh_in;
}


} // Discretization namespace
} // AMP namespace


