#ifndef included_MultiDOF_Manager
#define included_MultiDOF_Manager

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_Manager.h"
#include "discretization/DOF_ManagerParameters.h"


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

    using DOFManager::getDOFs;

    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new multiDOFManager object.
     * \param comm  Comm over which the DOFManager will exist
     * \param managers  List of the DOFManagers on the current processor
     */
    multiDOFManager ( AMP_MPI comm, std::vector<DOFManager::shared_ptr> managers );


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


    //! Get the DOFManagers that compose the multiDOFManager
    std::vector<DOFManager::shared_ptr>  getDOFManagers() const;


    /** \brief Subset the DOF Manager for a mesh
     * \details  This will subset a DOF manager for a particular mesh.  The resulting DOFManager
     *    can exist on either the comm of the parent DOF manager, or the comm of the mesh (default).
     * \param[in]  mesh         The mesh to use to subset
     * \param[in]  useMeshComm  Do we want to use the mesh comm for the new DOFManager.
     *                          Note: if this is true, any processors that do not contain the mesh will return NULL.
     */
    virtual DOFManager::shared_ptr subset( const AMP::Mesh::Mesh::shared_ptr mesh, bool useMeshComm=true );


    /** \brief Subset the DOF Manager for a mesh element iterator
     * \details  This will subset a DOF manager for a given mesh element iterator.  
     *    The resulting DOFManager will exist on the privided comm.
     * \param[in]  mesh     The mesh iterator for the subset
     * \param[in]  comm     The desired comm
     */
    virtual DOFManager::shared_ptr subset( const AMP::Mesh::MeshIterator &iterator, AMP_MPI comm );


private:
    std::vector<DOFManager::shared_ptr>                         d_managers;
    std::vector<size_t>                                         d_localSize;
    std::vector<size_t>                                         d_globalSize;

    // Data used to convert between the local (sub) and global (parent) DOFs
    struct subDOF_struct {
        size_t DOF1_begin;
        size_t DOF1_end;
        size_t DOF2_begin;
        size_t DOF2_end;
        // Constructors
        inline subDOF_struct( size_t v1, size_t v2, size_t v3, size_t v4 ) {
            DOF1_begin = v1;
            DOF1_end   = v2;
            DOF2_begin = v3;
            DOF2_end   = v4;
        } 
        inline subDOF_struct( ) {
            DOF1_begin = ~size_t(0);
            DOF1_end   = ~size_t(0);
            DOF2_begin = ~size_t(0);
            DOF2_end   = ~size_t(0);
        }
        // Overload key operators
        inline bool operator== (const subDOF_struct& rhs ) const {
            return DOF1_begin==rhs.DOF1_begin && DOF1_end==rhs.DOF1_end &&
                   DOF2_begin==rhs.DOF2_begin && DOF2_end==rhs.DOF2_end;
        }
        inline bool operator!= (const subDOF_struct& rhs ) const {
            return DOF1_begin!=rhs.DOF1_begin || DOF1_end!=rhs.DOF1_end ||
                   DOF2_begin!=rhs.DOF2_begin || DOF2_end!=rhs.DOF2_end;
        }
        inline bool operator>= (const subDOF_struct& rhs ) const {
            if ( DOF1_begin != rhs.DOF1_begin )
                return DOF1_begin>=rhs.DOF1_begin;
            if ( DOF1_end != rhs.DOF1_end )
                return DOF1_end>=rhs.DOF1_end;
            if ( DOF2_begin != rhs.DOF2_begin )
                return DOF2_begin>=rhs.DOF2_begin;
            return DOF2_end>=rhs.DOF2_end;
        }
        inline bool operator> (const subDOF_struct& rhs ) const {
            return operator>=(rhs) && operator!=(rhs);
        }
        inline bool operator< (const subDOF_struct& rhs ) const {
            return !operator>=(rhs);
        }
        inline bool operator<= (const subDOF_struct& rhs ) const {
            return !operator>=(rhs) || operator==(rhs);
        }
    };
    std::vector< std::vector<subDOF_struct> >      d_subToGlobalDOF;
    std::vector< std::vector<subDOF_struct> >      d_globalToSubDOF;
};


}
}

#endif

