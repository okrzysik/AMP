#ifndef included_simpleDOF_Manager
#define included_simpleDOF_Manager

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "discretization/DOF_ManagerParameters.h"


namespace AMP {
namespace Discretization {

/**
 * \class simpleDOFManager
 * \brief A derived class to create a simple DOF_Manager
 * \details  This derived class impliments a concrete DOF_Manager for creating Vectors 
 *    over a mesh on a particular mesh entity.  For example it can create a NodalVector
 *    over the entire Mesh.  Note: this class will be replaced by a more complete 
 *    Discretization interface.
 */
class simpleDOFManager: public DOFManager
{
public:

    using DOFManager::getDOFs;


    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new DOF manager object.
     * \param mesh      Mesh over which we want to construct the DOF map
     * \param type      The geometric entity type for the DOF map
     * \param gcw       The desired ghost width
     * \param split     Do we want to split the DOFManager by the meshes returning a multiDOFManager
     */
    static DOFManager::shared_ptr  create( boost::shared_ptr<AMP::Mesh::Mesh> mesh, 
        AMP::Mesh::GeomType type, int gcw, int DOFsPerElement, bool split=true );


    /**
     * \brief Create a new DOF manager object
     * \details  This is will create a new simpleDOFManager from a mesh iterator
     * \param mesh      Mesh over which the iterators are defined
     * \param it1       The iterator over the elements (including ghost cells)
     * \param it2       The iterator over the elements (excluding ghost cells)
     */
    static DOFManager::shared_ptr  create( boost::shared_ptr<AMP::Mesh::Mesh> mesh, 
        const AMP::Mesh::MeshIterator it1, const AMP::Mesh::MeshIterator it2, int DOFsPerElement );


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     *  Note: this function only works if the element we are search for is a element on which a DOF exists
     *  (the underlying mesh element type must match the geometric entity type specified at construction).
     * \param[in]  id       The element ID to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated
     *  with the DOFs.  Each element in the iterator will have 1 or more DOFs
     *  that are associated with that element.  For eaxample, a NodalVectorDOF
     *  would have 3 DOFs stored at each node, and would return an iterator over
     *  all the nodes. 
     */
    virtual AMP::Mesh::MeshIterator getIterator() const;


    //! Get the remote DOFs for a vector
    virtual std::vector<size_t> getRemoteDOFs() const;


    //! Get the row DOFs given a mesh element
    virtual std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const;


private:
    // Function to find the remote DOF given a set of mesh element IDs
    std::vector<size_t> getRemoteDOF(std::vector<AMP::Mesh::MeshElementID> remote_ids ) const;

    // Function to initialize the data
    void initialize();

    // Data members
    boost::shared_ptr<AMP::Mesh::Mesh>  d_mesh;
    AMP::Mesh::GeomType d_type;
    AMP::Mesh::MeshIterator d_localIterator;
    AMP::Mesh::MeshIterator d_ghostIterator;
    int DOFsPerElement;
    std::vector<AMP::Mesh::MeshElementID> d_local_id;
    std::vector<AMP::Mesh::MeshElementID> d_remote_id;
    std::vector<size_t> d_remote_dof;
};


}
}

#endif

