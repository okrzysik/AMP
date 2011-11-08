#ifndef included_AMP_DOF_ManagerParameters
#define included_AMP_DOF_ManagerParameters

#include <boost/shared_ptr.hpp>
#include "ampmesh/Mesh.h"

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

    typedef boost::shared_ptr<AMP::Discretization::DOFManagerParameters>  shared_ptr;

    //! Empty constructor for a DOF manager object
    DOFManagerParameters( );

    //! Default constructor for a DOF manager object
    DOFManagerParameters( boost::shared_ptr<AMP::Mesh::Mesh> mesh );

    //! Return the mesh
    boost::shared_ptr<AMP::Mesh::Mesh> getMesh() { return mesh; }


protected:

    //! Pointer to the underlying Mesh (may be NULL)
    boost::shared_ptr<AMP::Mesh::Mesh>  mesh;

    //! Pointer to the underlying VectorSpace (may be NULL)
    //boost::shared_ptr<AMP::Discretization::VectorSpace>  vectorSpace;

};



} // Discretization namespace
} // AMP namespace

#endif

