#ifndef included_AMP_CommVariable_H
#define included_AMP_CommVariable_H

#include "utils/AMP_MPI.h"
#include "vectors/SubsetVariable.h"

namespace AMP {
namespace LinearAlgebra {


/** \class MeshVariable
  * \brief An AMP Variable that describes how to subset a DOF for a mesh
  * \see SubsetVector
  */
class CommVariable : public SubsetVariable
{
public:
    /** \brief Constructor
      * \param[in] name  The name of the new variable
      * \param[in] comm  The AMP_MPI communicator of the new variable
      */
    CommVariable( const std::string &name, AMP_MPI comm );

    virtual AMP::Discretization::DOFManager::shared_ptr
        getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr ) const;

private:
    CommVariable();
    AMP_MPI d_comm;
};
}
}

#endif
