#ifndef included_AMP_DOF_Manager
#define included_AMP_DOF_Manager

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "discretization/DOF_ManagerParameters.h"
#include "utils/AMP_MPI.h"


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
class DOFManager: public boost::enable_shared_from_this<AMP::Discretization::DOFManager>
{
public:

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a DOF manager object.
     */
    typedef boost::shared_ptr<AMP::Discretization::DOFManager>  shared_ptr;

    /**
     *\typedef const_shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a DOF manager object.
     */
    typedef boost::shared_ptr<const AMP::Discretization::DOFManager>  const_shared_ptr;


    /** \brief Basic constructor for DOFManager
     * \details  This will create a very simple DOFManager with the given number of DOFs on each processor.
     *   It will not contain info to relate that to a mesh.  A derived implimentation should be used for
     *   more advanced features.  For example see simpleDOFManager and multiDOFManager.
     * \param[in]  N_local  The local number of DOFs
     * \param[in]  comm     The comm over which the DOFManager exists
     */
    DOFManager( size_t N_local, AMP_MPI comm );

    //! Deconstructor
    virtual ~DOFManager();

    /** \brief  Compares two DOFManager for equality.
      * \details This operation compares two DOF managers to see if they are equivalent
      * \param  rhs     DOFManager to compare 
      */
    virtual bool operator == ( const DOFManager & rhs ) const;


    /** \brief  Inverse of ==
      * \details This function performs an equality check and negates it.  Hence, it is not virtual
      * \param  rhs     DOFManager to compare
      */
    bool operator != ( const DOFManager & rhs ) const;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  id       The element ID to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <size_t> &dofs ) const;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  ids      The element IDs to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const std::vector<AMP::Mesh::MeshElementID> &ids, std::vector <size_t> &dofs ) const;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated
     *  with the DOFs.  Each element in the iterator will have 1 or more DOFs
     *  that are associated with that element.  For example, a NodalVectorDOF
     *  would have 3 DOFs stored at each node, and would return an iterator over
     *  all the nodes. 
     */
    virtual AMP::Mesh::MeshIterator getIterator() const;


    /** \brief  The first D.O.F. on this core
     * \return The first D.O.F. on this core
     */
    virtual size_t  beginDOF() const;


    /** \brief  One past the last D.O.F. on this core
     * \return One past the last D.O.F. on this core
     */
    virtual size_t  endDOF() const;


    /** \brief  The local number of D.O.F 
     * \return  The local number of D.O.F 
     */
    virtual size_t  numLocalDOF() const;


    /** \brief  The global number of D.O.F 
     * \return  The global number of D.O.F 
     */
    virtual size_t  numGlobalDOF() const;
 

    //! Get the comm for the DOFManger
    virtual AMP_MPI  getComm() const;
 

    //! Get the remote DOFs for a vector
    virtual std::vector<size_t> getRemoteDOFs() const;


    //! Get the row DOFs given a mesh element
    virtual std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const;


    /** \brief Subset the DOF Manager for a AMP_MPI communicator
     * \details  This will subset a DOF manager for a given communicator.
     * \param[in]  comm         The communicator to use to subset
     */
    virtual DOFManager::shared_ptr subset( AMP_MPI comm );


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
     * \param[in]  iterator     The mesh iterator for the subset
     * \param[in]  comm         The desired comm
     */
    virtual DOFManager::shared_ptr subset( const AMP::Mesh::MeshIterator &iterator, AMP_MPI comm );

protected:

    //!  Empty constructor for a DOF manager object
    DOFManager() {};

    //! The DOF manager parameters
    const DOFManagerParameters::shared_ptr params;

    //! The begining DOF, ending DOF and number of local DOFs for this processor
    size_t d_begin, d_end, d_global;

    //! The comm for this DOFManager
    AMP_MPI d_comm;
};



} // Discretization namespace
} // AMP namespace

#endif

