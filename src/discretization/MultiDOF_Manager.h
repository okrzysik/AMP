#ifndef included_AMP_MultiDOF_Manager
#define included_AMP_MultiDOF_Manager

#include "AMP/ampmesh/Mesh.h"
#include "AMP/ampmesh/MeshElement.h"
#include "AMP/discretization/DOF_Manager.h"
#include <memory>


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
class multiDOFManager : public DOFManager
{
public:
    using DOFManager::subset;

    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new multiDOFManager object.
     * \param comm  Comm over which the DOFManager will exist
     * \param managers  List of the DOFManagers on the current processor
     */
    multiDOFManager( const AMP_MPI &comm, std::vector<DOFManager::shared_ptr> managers );

    //! Deconstructor
    ~multiDOFManager() override;

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
    void getDOFs( const AMP::Mesh::MeshElementID &id, std::vector<size_t> &dofs ) const override;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  ids      The element IDs to collect nodal objects for.
     *                      Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    void getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids,
                  std::vector<size_t> &dofs ) const override;

    /** \brief Get the mesh element for a DOF
     * \details  This will return the mesh element associated with a given DOF.
     * \param[in] dof       The entry in the vector associated with DOF
     * @return              The element for the given DOF.
     */
    AMP::Mesh::MeshElement getElement( size_t dof ) const override;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated with the DOFs.
     * Note: if any sub-DOFManagers are the same, then this will iterate over repeated elements.
     */
    AMP::Mesh::MeshIterator getIterator() const override;


    //! Get the remote DOFs for a vector
    std::vector<size_t> getRemoteDOFs() const override;


    //! Get the row DOFs given a mesh element
    std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const override;


    /** \brief Subset the DOF Manager for a AMP_MPI communicator
     * \details  This will subset a DOF manager for a given communicator.
     * \param[in]  comm         The communicator to use to subset
     */
    DOFManager::shared_ptr subset( const AMP_MPI &comm ) override;


    /** \brief Subset the DOF Manager for a mesh
     * \details  This will subset a DOF manager for a particular mesh.  The resulting DOFManager
     *    can exist on either the comm of the parent DOF manager, or the comm of the mesh (default).
     * \param[in]  mesh         The mesh to use to subset
     * \param[in]  useMeshComm  Do we want to use the mesh comm for the new DOFManager.
     *                          Note: if this is true, any processors that do not contain the mesh
     * will return NULL.
     */
    DOFManager::shared_ptr subset( const AMP::Mesh::Mesh::shared_ptr mesh,
                                   bool useMeshComm = true ) override;


    /** \brief Subset the DOF Manager for a mesh element iterator
     * \details  This will subset a DOF manager for a given mesh element iterator.
     *    The resulting DOFManager will exist on the privided comm.
     * \param[in]  iterator     The mesh iterator for the subset
     * \param[in]  comm         The desired comm
     */
    DOFManager::shared_ptr subset( const AMP::Mesh::MeshIterator &iterator,
                                   const AMP_MPI &comm ) override;

public:
    //! Get the DOFManagers that compose the multiDOFManager
    std::vector<DOFManager::shared_ptr> getDOFManagers() const;


    /** \brief   Function to convert DOFs from a sub-manager DOF to the global DOF
     * \details  This function returns the global DOF given the local DOF.  Note that
     *      subDOFManager is specified by the index in the std::vector of DOFManagers
     *      (see getDOFManagers).  This is needed since the same DOFManager may be
     *      repeated many times so searching is not an option.  For example, consider
     *      a multiVector with multiple vectors of the same type.
     * \param[in]  DOFManager       The index to the desired DOFManager (see getDOFManagers)
     * \param[in]  localDOF         The local DOF to convert to global DOF
     */
    std::vector<size_t> getGlobalDOF( const int DOFManager,
                                      const std::vector<size_t> &localDOF ) const;


    /** Function to convert DOFs from the global DOF to a sub-manager DOF
     *  If a given global DOF is not in the given sub-manager, then -1 will
     *  be returned for its value.
     */
    /** \brief   Function to convert DOFs from the global DOF to a sub-manager DOF
     * \details  This function returns the DOF on a a sub-manager given the global DOF.
     *      If a given global DOF is not in the given sub-manager, then -1 will
     *      be returned for its value.  Note that subDOFManager is specified by the
     *      index in the std::vector of DOFManagers (see getDOFManagers).
     *      This is needed since the same DOFManager may be repeated many times so
     *      searching is not an option.  For example, consider a multiVector with
     *      multiple vectors of the same type.
     * \param[in]  DOFManager       The index to the desired DOFManager (see getDOFManagers)
     * \param[in]  globalDOF        The global DOF to convert to a DOF on the sub DOFManager
     */
    std::vector<size_t> getSubDOF( const int DOFManager, std::vector<size_t> &globalDOF ) const;


private:
    // Convert the local to global dof
    inline size_t subToGlobal( int manager, size_t dof ) const;

    // Convert the global to local dof
    inline std::pair<size_t, int> globalToSub( size_t dof ) const;

private:
    // Data used to convert between the local (sub) and global (parent) DOFs
    struct DOFMapStruct {
        // Constructors
        inline DOFMapStruct( size_t sub_start, size_t sub_end, size_t global_start, size_t id )
        {
            data[0] = sub_start;
            data[1] = sub_end;
            data[2] = global_start;
            data[3] = id;
        }
        inline DOFMapStruct()
        {
            data[0] = 0;
            data[1] = 0;
            data[2] = 0;
            data[3] = 0;
        }
        // Convert ids
        inline size_t toGlobal( size_t local ) const { return local - data[0] + data[2]; }
        inline size_t toLocal( size_t global ) const { return global - data[2] + data[0]; }
        inline bool inRangeLocal( size_t local ) const
        {
            return local >= data[0] && local < data[1];
        }
        inline size_t inRangeGlobal( size_t global ) const
        {
            return global >= data[2] && ( global - data[2] ) < ( data[1] - data[0] );
        }
        inline size_t id() const { return data[3]; }

    private:
        size_t data[4];
    };


private:
    std::vector<DOFManager::shared_ptr> d_managers;
    std::vector<size_t> d_ids;
    std::vector<size_t> d_localSize;
    std::vector<size_t> d_globalSize;
    std::vector<DOFMapStruct> d_dofMap;
    const size_t neg_one = ~( (size_t) 0 );
};
} // namespace Discretization
} // namespace AMP

#endif
