#ifndef included_AMP_AMPCSRMatrixParameters
#define included_AMP_AMPCSRMatrixParameters

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/MatrixParameters.h"
#include "AMP/utils/AMP_MPI.h"
#include "AMP/vectors/Vector.h"

namespace AMP::LinearAlgebra {

/** \class MatrixParameters
 * \brief  A class used to hold basic parameters for a matrix
 */
template<typename Config>
class AMPCSRMatrixParameters : public MatrixParameters
{
    using gidx_t = typename Config::gidx_t;
    using lidx_t = typename Config::lidx_t;

public:
    AMPCSRMatrixParameters() = delete;

    /** \brief Constructor
     * \param[in] dofLeft    The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left vector )
     * \param[in] dofRight   The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm       Communicator for the matrix
     * \param[in] getRowNNZ  Function that finds NNZ counts of each row
     * \param[in] getRowCols Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft    The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left vector )
     * \param[in] dofRight   The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm       Communicator for the matrix
     * \param[in] backend    Acceleration backend for matrix operations
     * \param[in] getRowNNZ  Function that finds NNZ counts of each row
     * \param[in] getRowCols Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        AMP::Utilities::Backend backend,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft    The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left vector )
     * \param[in] dofRight   The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm       Communicator for the matrix
     * \param[in] varLeft    Pointer to left variable
     * \param[in] varRight   Pointer to right variable
     * \param[in] getRowNNZ  Function that finds NNZ counts of each row
     * \param[in] getRowCols Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<Variable> varLeft,
        std::shared_ptr<Variable> varRight,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft    The DOFManager for the left vector ( For \f$\mathbf{y}^T\mathbf{Ax}\f$,
     * \f$y\f$ is a left vector )
     * \param[in] dofRight   The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm       Communicator for the matrix
     * \param[in] varLeft    Pointer to left variable
     * \param[in] varRight   Pointer to right variable
     * \param[in] backend    Acceleration backend for matrix operations
     * \param[in] getRowNNZ  Function that finds NNZ counts of each row
     * \param[in] getRowCols Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<Variable> varLeft,
        std::shared_ptr<Variable> varRight,
        AMP::Utilities::Backend backend,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft       The DOFManager for the left vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
     * \param[in] dofRight      The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm          Communicator for the matrix
     * \param[in] commListLeft  Communication list for the left vector
     * \param[in] commListRight Communication list for the right vector
     * \param[in] getRowNNZ     Function that finds NNZ counts of each row
     * \param[in] getRowCols    Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<CommunicationList> commListLeft,
        std::shared_ptr<CommunicationList> commListRight,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft       The DOFManager for the left vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
     * \param[in] dofRight      The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm          Communicator for the matrix
     * \param[in] commListLeft  Communication list for the left vector
     * \param[in] commListRight Communication list for the right vector
     * \param[in] backend       Acceleration backend for matrix operations
     * \param[in] getRowNNZ     Function that finds NNZ counts of each row
     * \param[in] getRowCols    Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<CommunicationList> commListLeft,
        std::shared_ptr<CommunicationList> commListRight,
        AMP::Utilities::Backend backend,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft       The DOFManager for the left vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
     * \param[in] dofRight      The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm          Communicator for the matrix
     * \param[in] varLeft       Pointer to left variable
     * \param[in] varRight      Pointer to right variable
     * \param[in] commListLeft  Communication list for the left vector
     * \param[in] commListRight Communication list for the right vector
     * \param[in] getRowNNZ     Function that finds NNZ counts of each row
     * \param[in] getRowCols    Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<Variable> varLeft,
        std::shared_ptr<Variable> varRight,
        std::shared_ptr<CommunicationList> commListLeft,
        std::shared_ptr<CommunicationList> commListRight,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    /** \brief Constructor
     * \param[in] dofLeft       The DOFManager for the left vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$y\f$ is a left vector )
     * \param[in] dofRight      The DOFManager for the right vector ( For
     * \f$\mathbf{y}^T\mathbf{Ax}\f$, \f$x\f$ is a right vector )
     * \param[in] comm          Communicator for the matrix
     * \param[in] varLeft       Pointer to left variable
     * \param[in] varRight      Pointer to right variable
     * \param[in] commListLeft  Communication list for the left vector
     * \param[in] commListRight Communication list for the right vector
     * \param[in] backend       Acceleration backend for matrix operations
     * \param[in] getRowNNZ     Function that finds NNZ counts of each row
     * \param[in] getRowCols    Function that finds finds column indices of a row in each block
     */
    explicit AMPCSRMatrixParameters(
        std::shared_ptr<AMP::Discretization::DOFManager> dofLeft,
        std::shared_ptr<AMP::Discretization::DOFManager> dofRight,
        const AMP_MPI &comm,
        std::shared_ptr<Variable> varLeft,
        std::shared_ptr<Variable> varRight,
        std::shared_ptr<CommunicationList> commListLeft,
        std::shared_ptr<CommunicationList> commListRight,
        AMP::Utilities::Backend backend,
        const std::function<void( const gidx_t, lidx_t &, lidx_t & )> getRowNNZ  = {},
        const std::function<void( const gidx_t, gidx_t *, gidx_t * )> getRowCols = {} );

    //! Get the bound function that finds NNZ counts each row
    const std::function<void( const gidx_t, lidx_t &, lidx_t & )> &getRowNNZFunction() const
    {
        return d_getRowNNZ;
    }

    //! Get the bound function that finds column indices of row in each block
    const std::function<void( const gidx_t, gidx_t *, gidx_t * )> &getRowColsFunction() const
    {
        return d_getRowCols;
    }

protected:
    //! Function to get number of non-zeros per row in each block
    std::function<void( const gidx_t, lidx_t &, lidx_t & )> d_getRowNNZ;

    //! Function to get column indices of row in each block
    std::function<void( const gidx_t, gidx_t *, gidx_t * )> d_getRowCols;
};

} // namespace AMP::LinearAlgebra

#endif
