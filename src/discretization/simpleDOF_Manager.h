#ifndef included_simpleDOF_Manager
#define included_simpleDOF_Manager

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
 * \class simpleDOF_Manager
 * \brief A derived class to create a simple DOF_Manager
 * \details  This derived class impliments a concrete DOF_Manager for creating Vectors 
 *    over a mesh on a particular mesh entity.  For example it can create a NodalVector
 *    over the entire Mesh.  Note: this class will be replaced by a more complete 
 *    Discretization interface.
 */
class simpleDOFManager: public DOFManager, public boost::enable_shared_from_this<AMP::Discretization::simpleDOFManager>
{
public:

    /**
     * \brief Create a new DOF manager object
     * \details  This is the standard constructor for creating a new DOF manager object.
     * \param mesh  Mesh over which we want to construct the DOF map
     * \param type  The geometric entity type for the DOF map
     * \param gcw   The desired ghost width
     */
    simpleDOFManager ( boost::shared_ptr<AMP::Mesh::Mesh> mesh, AMP::Mesh::GeomType type, int gcw, int DOFsPerElement );


    /** \brief Get the entry indices of DOFs given a mesh element
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     * \param[in]  obj      The element to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     * \param[in]  which    Which D.O.F. to get.  If not specified, return all D.O.F.s
     */
    virtual void getDOFs( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &dofs , std::vector<unsigned int> which = std::vector<unsigned int>(0) ) const;


    /** \brief Get the entry indices of DOFs given a mesh element ID
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     *  Note: this function only works if the element we are search for is a element on which a DOF exists
     *  (the underlying mesh element type must match the geometric entity type specified at construction).
     * \param[in]  id       The element ID to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] dofs     The entries in the vector associated with D.O.F.s on the nodes
     */
    virtual void getDOFs( const AMP::Mesh::MeshElementID &id, std::vector <unsigned int> &dofs ) const;


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
    virtual size_t  beginDOF ( ) const;


    /** \brief  One past the last D.O.F. on this core
     * \return One past the last D.O.F. on this core
     */
    virtual size_t  endDOF ( ) const;


    /** \brief  The local number of D.O.F 
     * \return  The local number of D.O.F 
     */
    virtual size_t  numLocalDOF ( ) const;


    /** \brief  The global number of D.O.F 
     * \return  The global number of D.O.F 
     */
    virtual size_t  numGlobalDOF ( ) const;
 

    //! Get the comm for the DOFManger
    virtual AMP_MPI  getComm() const;
 

    //! Get the remote DOFs for a vector
    virtual std::vector<size_t> getRemoteDOFs() const;


    //! Get the row DOFs given a mesh element
    virtual std::vector<size_t> getRowDOFs( const AMP::Mesh::MeshElement &obj ) const;


    /**
     * \brief Create a new AMP vector
     * \details  This function creates a new AMP vector for the given variable, using the current DOF properties.
     * \param variable  Variable that will be used to create the vector
     */
    AMP::LinearAlgebra::Vector::shared_ptr   createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable );


    /**
     * \brief Create a new AMP matrix
     * \details  This function creates a new AMP matrix for the given variable, using the current DOF properties.
     * \param operand  Variable that will be used to create the matrix
     * \param result   Variable that will be used to create the matrix
     */
    AMP::LinearAlgebra::Matrix::shared_ptr   createMatrix ( AMP::LinearAlgebra::Vector::shared_ptr operand , AMP::LinearAlgebra::Vector::shared_ptr result );


private:
    // Function to find the remote DOF given a set of mesh element IDs
    std::vector<size_t> getRemoteDOF(std::vector<AMP::Mesh::MeshElementID> remote_ids ) const;

    // Data members
    boost::shared_ptr<AMP::Mesh::Mesh>  d_mesh;
    AMP::Mesh::GeomType d_type;
    int d_gcw;
    int DOFsPerElement;
    std::vector<AMP::Mesh::MeshElementID> d_local_id;
    std::vector<AMP::Mesh::MeshElementID> d_remote_id;
    std::vector<size_t> d_remote_dof;
    size_t d_begin, d_end, d_global;
};


}
}

#endif

