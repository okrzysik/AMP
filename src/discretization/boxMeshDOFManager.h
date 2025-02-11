#ifndef included_AMP_boxMeshDOFManager
#define included_AMP_boxMeshDOFManager

#include "AMP/discretization/simpleDOF_Manager.h"
#include "AMP/mesh/structured/BoxMesh.h"


namespace AMP::Discretization {


/**
 * \class boxMeshDOFManager
 * \brief A derived class to create a simple DOF manager on a box mesh
 * \details  This derived class implements a simpleDOFManager for creating Vectors
 *    over a mesh on boxMesh.
 */
class boxMeshDOFManager : public simpleDOFManager
{
public:
    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new DOF manager object.
     * \param mesh          Mesh over which we want to construct the DOF map (must be a boxMesh)
     * \param type          The geometric entity type for the DOF map
     * \param gcw           The desired ghost width
     * \param DOFsPerElement The desired number of DOFs per element
     * multiDOFManager
     */
    boxMeshDOFManager( std::shared_ptr<const AMP::Mesh::Mesh> mesh,
                       AMP::Mesh::GeomType type,
                       int gcw,
                       int DOFsPerElement );


    //! Destructor
    virtual ~boxMeshDOFManager() = default;

    //! Get the underlying box mesh
    auto getBoxMesh() const { return d_boxMesh; }

    //! Get the array size for the local variables
    ArraySize getArraySize() const;

private:
    // Private constructor
    boxMeshDOFManager() = delete;

private: // Data
    std::shared_ptr<const AMP::Mesh::BoxMesh> d_boxMesh;
};

} // namespace AMP::Discretization

#endif
