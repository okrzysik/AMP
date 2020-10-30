#ifndef included_AMP_ManagedEpetraMatrix
#define included_AMP_ManagedEpetraMatrix

#include <set>

#include "AMP/discretization/DOF_Manager.h"
#include "AMP/matrices/ManagedMatrix.h"
#include "AMP/matrices/trilinos/EpetraMatrix.h"
#include "AMP/matrices/trilinos/ManagedEpetraMatrixParameters.h"
#include "AMP/vectors/trilinos/epetra/EpetraVector.h"

#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {


/** \class ManagedEpetraMatrix
 * \brief  A class that wraps an Epetra_CrsMatrix
 * \details  This class stores an Epetra_FECrsMatrix and provides
 * the AMP interface to this matrix.
 */
class ManagedEpetraMatrix : public EpetraMatrix, public ManagedMatrix
{
protected:
    //!  Parameters used to construct the matrix
    std::shared_ptr<ManagedEpetraMatrixParameters> d_pParameters;

    //!  Empty constructor
    ManagedEpetraMatrix() = delete;

    //!  Copy constructor
    ManagedEpetraMatrix( const ManagedEpetraMatrix &rhs );

    //!  Assignment operator
    ManagedEpetraMatrix &operator=( const ManagedEpetraMatrix &rhs ) = delete;

    //!  \f$A_{i,j}\f$ storage of off-core data
    std::map<int, std::map<int, double>> d_OtherData;

    //!  Update data off-core
    void setOtherData();

    void multiply( shared_ptr other_op, shared_ptr &result ) override;

public:
    /** \brief Constructor
     * \param[in] p  The description of the matrix
     */
    explicit ManagedEpetraMatrix( std::shared_ptr<ManagedEpetraMatrixParameters> p );

    /** \brief Constructor from Epetra_CrsMatrix
     * \param[in]  m  Matrix to wrap
     * \param[in]  dele  If true, this class deletes the matrix
     */
    explicit ManagedEpetraMatrix( Epetra_CrsMatrix *m, bool dele = false );

    //! Destructor
    virtual ~ManagedEpetraMatrix() {}

    void createValuesByGlobalID( size_t, const std::vector<size_t> & ) override;


    void mult( const Vector::const_shared_ptr in, Vector::shared_ptr out ) override;
    void multTranspose( const Vector::const_shared_ptr in, Vector::shared_ptr out ) override;
    Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const override;
    void scale( double alpha ) override;
    void axpy( double alpha, const Matrix &rhs ) override;
    size_t numGlobalRows() const override { return d_epetraMatrix->NumGlobalRows(); }
    size_t numGlobalColumns() const override { return d_epetraMatrix->NumGlobalCols(); }
    void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    void getRowByGlobalID( size_t row,
                           std::vector<size_t> &cols,
                           std::vector<double> &values ) const override;


    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    std::vector<size_t> getColumnIDs( size_t row ) const override;

    void getValuesByGlobalID( size_t num_rows,
                              size_t num_cols,
                              size_t *rows,
                              size_t *cols,
                              double *values ) const override;

    void setScalar( double ) override;
    void setDiagonal( Vector::const_shared_ptr in ) override;

    void makeConsistent() override;
    double L1Norm() const override;
    Matrix::shared_ptr cloneMatrix() const override;
    Vector::shared_ptr getRightVector() const override;
    Vector::shared_ptr getLeftVector() const override;
    Discretization::DOFManager::shared_ptr getRightDOFManager() const override;
    Discretization::DOFManager::shared_ptr getLeftDOFManager() const override;
    void fillComplete() override;
    void setIdentity() override;
    void zero() override;
};
} // namespace LinearAlgebra
} // namespace AMP

#include "ManagedEpetraMatrix.inline.h"
#endif
