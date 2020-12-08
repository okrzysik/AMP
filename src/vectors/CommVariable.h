#ifndef included_AMP_CommVariable_H
#define included_AMP_CommVariable_H

#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/SubsetVariable.h"

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
    CommVariable( const std::string &name, const AMP_MPI &comm );

    virtual AMP::Discretization::DOFManager::shared_ptr
        getSubsetDOF( std::shared_ptr<AMP::Discretization::DOFManager> ) const override;

private:
    CommVariable();
    AMP_MPI d_comm;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
