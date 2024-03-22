#ifndef included_AMP_MatrixParametersBase
#define included_AMP_MatrixParametersBase

#include "AMP/utils/AMP_MPI.h"
#include "AMP/utils/Utilities.h"
#include "AMP/vectors/Variable.h"

namespace AMP::LinearAlgebra {


/** \class MatrixParametersBase
 * \brief  A class used to hold basic parameters for a matrix
 */
class MatrixParametersBase
{
public:
    MatrixParametersBase() = delete;

    /** \brief Constructor
     * \param[in] comm     Communicator for the matrix
     */
    explicit MatrixParametersBase( const AMP_MPI &comm ) : d_comm( comm ) {}

    //! Deconstructor
    virtual ~MatrixParametersBase() = default;

    //!  Get the communicator for the matrix
    AMP::AMP_MPI &getComm() { return d_comm; }

    //! memory space where the matrix should live
    AMP::Utilities::MemoryType d_memory_location = AMP::Utilities::MemoryType::host;

    //!  The variable for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left
    //!  vector )
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableLeft;

    //!  The variable for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right
    //!  vector )
    std::shared_ptr<AMP::LinearAlgebra::Variable> d_VariableRight;

protected:
    // The comm of the matrix
    AMP_MPI d_comm;
};
} // namespace AMP::LinearAlgebra

#endif
