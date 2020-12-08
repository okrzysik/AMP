#ifndef included_AMP_MeshVariable_H
#define included_AMP_MeshVariable_H
#ifdef USE_AMP_MESH

#include "SubsetVariable.h"

namespace AMP {
namespace LinearAlgebra {


/** \class MeshVariable
 * \brief An AMP Variable that describes how to subset a DOF for a mesh
 * \see SubsetVector
 */
class MeshVariable : public SubsetVariable
{
public:
    /** \brief Constructor
     * \param[in] name  The name of the new variable
     * \param[in] mesh  The mesh to subset for
     * \param[in] useMeshComm  Use the comm of the mesh (otherwise use the comm of the parent
     * DOFManager)
     */
    MeshVariable( const std::string &name,
                  AMP::Mesh::Mesh::shared_ptr mesh,
                  bool useMeshComm = true );

    virtual AMP::Discretization::DOFManager::shared_ptr
        getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> ) const override;

private:
    MeshVariable();
    bool d_useMeshComm;
    AMP::Mesh::Mesh::shared_ptr d_mesh;
};


/** \class MeshItertorVariable
 * \brief An AMP Variable that describes how to subset a DOF for a mesh iterator
 * \see SubsetVector
 */
class MeshIteratorVariable : public SubsetVariable
{
public:
    /** \brief Constructor
     * \param[in] name     The name of the new variable
     * \param[in] iterator The iterator over the mesh elements of interest
     * \param[in] comm     The communicator to subset for
     */
    MeshIteratorVariable( const std::string &name,
                          const AMP::Mesh::MeshIterator &iterator,
                          const AMP_MPI &comm );

    virtual AMP::Discretization::DOFManager::shared_ptr
        getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> ) const override;

private:
    MeshIteratorVariable();
    const AMP::AMP_MPI d_comm;
    const AMP::Mesh::MeshIterator d_iterator;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
#endif
