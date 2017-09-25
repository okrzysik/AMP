#ifndef included_AMP_StridedVariable_H
#define included_AMP_StridedVariable_H

#include "SubsetVariable.h"

namespace AMP {
namespace LinearAlgebra {

/** \class StridedVariable
 * \brief An AMP Variable that describes how to stride a vector to create
 * a SubsetVector
 * \see SubsetVector
 * \see StridedIndexer
 */
class StridedVariable : public SubsetVariable
{
public:
    /** \brief Constructor
     * \param[in] name  The name of the new variable
     * \param[in] offset  The offset to start striding a vector
     * \param[in] stride  The stride of the vector
     */
    StridedVariable( const std::string &name, size_t offset, size_t stride );

    virtual AMP::Discretization::DOFManager::shared_ptr
        getSubsetDOF( AMP::Discretization::DOFManager::shared_ptr ) const override;

private:
    StridedVariable();
    size_t d_offset;
    size_t d_stride;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
