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
class DOFManager
{
public:

    /**
     *\typedef shared_ptr
     *\brief  Name for the shared pointer.
     *\details  Use this typedef for a reference counted pointer to a DOF manager object.
     */
    typedef boost::shared_ptr<AMP::Discretization::DOFManager>  shared_ptr;


    /** \brief Get the entry indices of DOFs given a mesh element
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  obj      The element to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] ids      The entries in the vector associated with D.O.F.s on the nodes
     * \param[in]  which    Which D.O.F. to get.  If not specified, return all D.O.F.s
     */
    virtual void getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids , std::vector<unsigned int> which = std::vector<unsigned int>(0) ) const;


    /** \brief   Get an entry over the mesh elements associated with the DOFs
     * \details  This will return an iterator over the mesh elements associated
     *  with the DOFs.  Each element in the iterator will have 1 or more DOFs
     *  that are associated with that element.  For eaxample, a NodalVectorDOF
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


protected:

    //!  Empty constructor for a DOF manager object
    DOFManager() {};

    //! The DOF manager parameters
    const DOFManagerParameters::shared_ptr params;

};



} // Discretization namespace
} // AMP namespace

#endif

