#ifndef included_MultiDOF_Manager
#define included_MultiDOF_Manager

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "discretization/DOF_ManagerParameters.h"
#include "vectors/Vector.h"
#include "matrices/Matrix.h"


namespace AMP {
namespace Discretization {


/**
 * \class multiDOFManager
 * \brief A derived class to combine multiple DOFManagers
 * \details  This derived class impliments a concrete DOF_Manager for creating DOFs that
 *   consist of multiple DOFManagers.  This is useful to combine multiple DOFs over meshes
 *   on a multiVector, for combining multiple discretizations, and for combining vectors.
 *   A multivector will have a pointer to a multiDOFManager instead of a standard DOFManager.
 *   It is also possible that a standard vector can use a multiDOFManager.
 */
class multiDOFManager: public DOFManager
{
public:

    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new multiDOFManager object.
     * \param comm  Comm over which the DOFManager will exist
     * \param managers  List of the DOFManagers on the current processor
     */
    multiDOFManager ( AMP_MPI comm, std::vector<DOFManager::shared_ptr> managers );


    /** \brief Get the entry indices of DOFs given a mesh element
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  obj      The element to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     * \param[in]  which    Which D.O.F. to get.  If not specified, return all D.O.F.s
     */
    virtual void getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <size_t> &dofs , std::vector<size_t> which = std::vector<size_t>(0) ) const;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     *  Note: this function only works if the element we are search for is a element on which a DOF exists
     *  (the underlying mesh element type must match the geometric entity type specified at construction).
     * \param[in]  id       The element ID to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated with the DOFs.  
     * Note: if any sub-DOFManagers are the same, then this will iterate over repeated elements.
     */
    virtual AMP::Mesh::MeshIterator getIterator() const;
 

    //! Get the remote DOFs for a vector
    virtual std::vector<size_t> getRemoteDOFs() const;


    //! Get the row DOFs given a mesh element
    virtual std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const;

    //! Function to convert DOFs from a sub-manager DOF to the global DOF
    std::vector<size_t>  getGlobalDOF(DOFManager::shared_ptr, std::vector<size_t>&) const;

    /** Function to convert DOFs from the global DOF to a sub-manager DOF
     *  If a given global DOF is not in the given sub-manager, then -1 will
     *  be returned for its value.
     */
    std::vector<size_t>  getSubDOF(DOFManager::shared_ptr, std::vector<size_t>&) const;

private:
    std::vector<DOFManager::shared_ptr>                         d_managers;
    std::vector<size_t>                                         d_localSize;
    std::vector<size_t>                                         d_globalSize;
    std::vector< std::vector< std::pair<size_t,size_t> > >      d_subToGlobalDOF;
    std::vector< std::vector< std::pair<size_t,size_t> > >      d_globalToSubDOF;
};


}
}

#endif

