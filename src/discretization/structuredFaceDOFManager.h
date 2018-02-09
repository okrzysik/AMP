#ifndef included_structuredFaceDOFManager
#define included_structuredFaceDOFManager

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/discretization/DOF_Manager.h"
#include "AMP/utils/shared_ptr.h"


namespace AMP {
namespace Discretization {


/**
 * \class structuredFaceDOFManager
 * \brief A derived class to create a DOFManager for faces
 * \details  This derived class impliments a concrete DOFManager for creating Vectors
 *    and matricies over a mesh on the faces of structured meshes.
 *    This is a specific implimentation designed for rectangular 3d meshes,
 *    and will create the unknowns on the faces.  Two faces are neighbors if they
 *    share an element.
 */
class structuredFaceDOFManager final : public DOFManager
{
public:
    using DOFManager::subset;


    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new DOF manager object.
     * \param mesh          Mesh over which we want to construct the DOF map
     * \param DOFsPerFace   The desired number of DOFs per face (x,y,z)
     * \param gcw           The desired ghost width (based on the volumes)
     */
    static DOFManager::shared_ptr
    create( AMP::shared_ptr<AMP::Mesh::Mesh> mesh, int DOFsPerFace[3], int gcw );


    //! Deconstructor
    virtual ~structuredFaceDOFManager();


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     *  Note: this function only works if the element we are search for is a element on which a DOF
     * exists
     *  (the underlying mesh element type must match the geometric entity type specified at
     * construction).
     * \param[in]  id       The element ID to collect nodal objects for.  Note: the mesh element may
     * be any type
     * (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const AMP::Mesh::MeshElementID &id,
                          std::vector<size_t> &dofs ) const override;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  ids      The element IDs to collect nodal objects for.
     *                      Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids,
                          std::vector<size_t> &dofs ) const override;


    /** \brief Get the mesh element for a DOF
     * \details  This will return the mesh element associated with a given DOF.
     * \param[in] dof       The entry in the vector associated with DOF
     * @return              The element for the given DOF.
     */
    virtual AMP::Mesh::MeshElement getElement( size_t dof ) const override;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated
     *  with the DOFs.  Each element in the iterator will have 1 or more DOFs
     *  that are associated with that element.  For eaxample, a NodalVectorDOF
     *  would have 3 DOFs stored at each node, and would return an iterator over
     *  all the nodes.
     */
    virtual AMP::Mesh::MeshIterator getIterator() const override;


    //! Get the remote DOFs for a vector
    virtual std::vector<size_t> getRemoteDOFs() const override;


    //! Get the row DOFs given a mesh element
    virtual std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const override;


private:
    // Empty constructor
    structuredFaceDOFManager() {}


private:
    // Function to find the remote DOF given a set of mesh element IDs
    std::vector<size_t>
    getRemoteDOF( const std::vector<AMP::Mesh::MeshElementID> &remote_ids ) const;

    // Function to initialize the data
    void initialize();

    // Append DOFs
    inline void appendDOFs( const AMP::Mesh::MeshElementID &id, std::vector<size_t> &dofs ) const;

    // Data members
    AMP::shared_ptr<AMP::Mesh::Mesh> d_mesh;
    int d_DOFsPerFace[3];
    int d_gcw;

    std::vector<AMP::Mesh::MeshElementID> d_local_ids[3];
    std::vector<AMP::Mesh::MeshElementID> d_remote_ids[3];
    std::vector<size_t> d_local_dofs[3];
    std::vector<size_t> d_remote_dofs[3];
};
} // namespace Discretization
} // namespace AMP

#endif
