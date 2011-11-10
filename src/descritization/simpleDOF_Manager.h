#ifndef included_simpleDOF_Manager
#define included_simpleDOF_Manager

#include "ampmesh/Mesh.h"
#include "ampmesh/MeshElement.h"
#include "descritization/DOF_Manager.h"
#include "descritization/DOF_ManagerParameters.h"

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
class simpleDOFManager: public DOFManager
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


    /** \brief Get the entry indices of nodal values given a mesh element
     * \param[in]  obj  The element to collect nodal objects for.  Note: the mesh element may be any type (include a vertex).
     * \param[out] ids  The entries in the vector associated with D.O.F.s on the nodes
     * \param[in]  which  Which D.O.F. to get.  If not specified, return all D.O.F.s
     * \details  This will return a vector of pointers into a Vector that are associated with which.
     */
    virtual void getDOFs ( const AMP::Mesh::MeshElement &obj, std::vector <unsigned int> &ids , unsigned int which = static_cast<unsigned int>(-1) ) const;


    /** \brief  The first D.O.F. on this core
     * \return The first D.O.F. on this core
     */
    virtual size_t  beginDOF ( );


    /** \brief  One past the last D.O.F. on this core
     * \return One past the last D.O.F. on this core
     */
    virtual size_t  endDOF ( );


    /**
     * \brief Create a new AMP vector
     * \details  This function creates a new AMP vector for the given variable, using the current DOF properties.
     * \param variable  Variable that will be used to create the vector
     */
    virtual AMP::LinearAlgebra::Vector::shared_ptr   createVector ( AMP::LinearAlgebra::Variable::shared_ptr variable );


    /**
     * \brief Create a new AMP matrix
     * \details  This function creates a new AMP matrix for the given variable, using the current DOF properties.
     * \param operand  Variable that will be used to create the matrix
     * \param result   Variable that will be used to create the matrix
     */
    //virtual  AMP::LinearAlgebra::Matrix::shared_ptr   createMatrix ( AMP::LinearAlgebra::Variable::shared_ptr operand , AMP::LinearAlgebra::Variable::shared_ptr result = AMP::LinearAlgebra::Variable::shared_ptr() );


private:
    boost::shared_ptr<AMP::Mesh::Mesh>  d_mesh;
    AMP::Mesh::GeomType d_type;
    int d_gcw;
    int DOFsPerElement;
    std::vector<AMP::Mesh::MeshElementID> d_local_id;
    std::vector<AMP::Mesh::MeshElementID> d_remote_id;
    std::vector<size_t> d_remote_dof;
    size_t d_begin, d_end;
};


}
}

#endif

