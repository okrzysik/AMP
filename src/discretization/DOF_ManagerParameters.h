#ifndef included_AMP_DOF_ManagerParameters
#define included_AMP_DOF_ManagerParameters

#include "AMP/ampmesh/Mesh.h"
#include <memory>

namespace AMP {
namespace Discretization {


/**
 * \class DOF_Manager
 * \brief A class used to provide DOF and vector creation routines
 *
 * \details  This class provides routines for calculating, accessing, and
 *    using the degrees of freedom (DOF) per object.  It is also responsible
 *    for creating vectors.
 */
class DOFManagerParameters
{
public:
    //! Empty constructor for a DOF manager object
    DOFManagerParameters();

    //! Default constructor for a DOF manager object
    explicit DOFManagerParameters( std::shared_ptr<AMP::Mesh::Mesh> mesh );

    //! Return the mesh
    std::shared_ptr<AMP::Mesh::Mesh> getMesh() { return mesh; }


protected:
    //! Pointer to the underlying Mesh (may be NULL)
    std::shared_ptr<AMP::Mesh::Mesh> mesh;

    //! Pointer to the underlying VectorSpace (may be NULL)
    // std::shared_ptr<AMP::Discretization::VectorSpace>  vectorSpace;
};


} // namespace Discretization
} // namespace AMP

#endif
