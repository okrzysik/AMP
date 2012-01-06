#ifndef included_subsetDOF_Manager
#define included_subsetDOF_Manager

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "discretization/DOF_ManagerParameters.h"


namespace AMP {
namespace Discretization {


/**
 * \class subsetDOFManager
 * \brief A derived class to subset a DOFManagers
 * \details  This derived class impliments a concrete DOF_Manager for maintaining
 *   a subset of a DOFManager.
 */
class subsetDOFManager: public DOFManager
{
public:

    using DOFManager::getDOFs;


    /** \brief Default constructor
     * \details  This is the default constructor for creating a subset DOF manager.  
     * \param[in]  parentDOFManager  The parent DOF manager
     * \param[out] dofs     The DOFs that will be part of the subset
     */
    subsetDOFManager( DOFManager::shared_ptr parentDOFManager, const std::vector <size_t> &dofs );


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


    //! Function to return the local DOFs on the parent DOF manager
    std::vector<size_t>  getLocalParentDOFs( ) const;


    //! Function to convert DOFs from a subset DOFManager DOF to the parent DOF
    std::vector<size_t>  getParentDOF( const std::vector<size_t>& ) const;


    /**
      *  Function to convert DOFs from the parent DOF to a subset manager DOF.
      *  Note: if the parent DOF does not exist in the subset, then -1 will be
      *  returned in it's place
      */
    std::vector<size_t>  getSubsetDOF( const std::vector<size_t>& ) const;


    //! Get the parent DOFManager
    DOFManager::shared_ptr  getDOFManager() const;


private:

    //! The parent DOF Manager
    DOFManager::shared_ptr d_parentDOFManager;

    //! The parent begin, end, and global DOFs
    size_t d_parentBegin, d_parentEnd, d_parentGlobal;

    //! The list of local DOFs (sorted, using the parent DOF numbering)
    std::vector<size_t> d_localDOFs;

    //! The list of remote DOFs
    std::vector<size_t> d_remoteParentDOFs;
    std::vector<size_t> d_remoteSubsetDOFs;

};


}
}

#endif

