#ifndef  included_AMP_MeshVariable_H
#define  included_AMP_MeshVariable_H
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
      */
    MeshVariable ( const std::string &name, AMP::Mesh::Mesh::shared_ptr mesh );

    virtual boost::shared_ptr<AMP::Discretization::subsetDOFManager>  getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr );

private:
    MeshVariable ();
    AMP::Mesh::Mesh::shared_ptr  d_mesh;
};


/** \class MeshItertorVariable
  * \brief An AMP Variable that describes how to subset a DOF for a mesh iterator
  * \see SubsetVector
  */
class MeshIteratorVariable : public SubsetVariable
{
public:
    /** \brief Constructor
      * \param[in] name  The name of the new variable
      * \param[in] mesh  The mesh to subset for
      */
    MeshIteratorVariable ( const std::string &name, const AMP::Mesh::MeshIterator &iterator );

    virtual boost::shared_ptr<AMP::Discretization::subsetDOFManager>  getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr );

private:
    MeshIteratorVariable ();
    const AMP::Mesh::MeshIterator  d_iterator;
};


}
}

#endif
#endif

