#ifndef included_AMP_MatrixParameters
#define included_AMP_MatrixParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/vectors/Vector.h"

namespace AMP {
namespace LinearAlgebra {


/** \class MatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
class MatrixParameters
{
public:
    //! Convenience typedef
    typedef AMP::shared_ptr<MatrixParameters> shared_ptr;

    /** \brief Constructor
     * \param[in] left     The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left
     * vector )
     * \param[in] right    The DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$x\f$ is a right
     * vector )
     * \param[in] comm     Communicator for the matrix
     */
    explicit MatrixParameters( AMP::Discretization::DOFManager::shared_ptr left,
                               AMP::Discretization::DOFManager::shared_ptr right,
                               const AMP_MPI &comm );

    //! Deconstructor
    virtual ~MatrixParameters(){};

    //! Return the local number of rows
    inline size_t getLocalNumberOfRows() const { return d_DOFManagerLeft->numLocalDOF(); }

    //! Return the local number of columns
    inline size_t getLocalNumberOfColumns() const { return d_DOFManagerRight->numLocalDOF(); }

    //! Return the global number of rows
    inline size_t getGlobalNumberOfRows() const { return d_DOFManagerLeft->numGlobalDOF(); }

    //! Return the global number of columns
    inline size_t getGlobalNumberOfColumns() const { return d_DOFManagerRight->numGlobalDOF(); }

    //!  Get the DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a
    //!  left vector )
    inline AMP::Discretization::DOFManager::shared_ptr getLeftDOFManager()
    {
        return d_DOFManagerLeft;
    }

    //!  Get the DOFManager for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a
    //!  right vector )
    inline AMP::Discretization::DOFManager::shared_ptr getRightDOFManager()
    {
        return d_DOFManagerRight;
    }

    //!  Get the communicator for the matric
    inline AMP::AMP_MPI &getComm() { return d_comm; }

    //!  The communication list of a left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a
    //!  left vector )
    CommunicationList::shared_ptr d_CommListLeft;

    //!  The communication list of a right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a
    //!  right vector )
    CommunicationList::shared_ptr d_CommListRight;

    //!  The variable for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left
    //!  vector )
    AMP::LinearAlgebra::Variable::shared_ptr d_VariableLeft;

    //!  The variable for the right vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right
    //!  vector )
    AMP::LinearAlgebra::Variable::shared_ptr d_VariableRight;

protected:
    MatrixParameters(){};

    // The DOFManager for the left vector ( may be null )
    AMP::Discretization::DOFManager::shared_ptr d_DOFManagerLeft;

    // The DOFManager for the right vector ( may be null )
    AMP::Discretization::DOFManager::shared_ptr d_DOFManagerRight;

    // The comm of the matrix
    AMP_MPI d_comm;
};
} // namespace LinearAlgebra
} // namespace AMP

#endif
