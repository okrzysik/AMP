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

    virtual void multiply( shared_ptr other_op, shared_ptr &result ) override;

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

    virtual void createValuesByGlobalID( size_t, const std::vector<size_t> & ) override;


    virtual void mult( const Vector::const_shared_ptr in, Vector::shared_ptr out ) override;
    virtual void multTranspose( const Vector::const_shared_ptr in,
                                Vector::shared_ptr out ) override;
    virtual Vector::shared_ptr
    extractDiagonal( Vector::shared_ptr buf = Vector::shared_ptr() ) const override;
    virtual void scale( double alpha ) override;
    virtual void axpy( double alpha, const Matrix &rhs ) override;
    virtual size_t numGlobalRows() const override { return d_epetraMatrix->NumGlobalRows(); }
    virtual size_t numGlobalColumns() const override { return d_epetraMatrix->NumGlobalCols(); }
    virtual void addValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    virtual void setValuesByGlobalID(
        size_t num_rows, size_t num_cols, size_t *rows, size_t *cols, double *values ) override;
    virtual void getRowByGlobalID( size_t row,
                                   std::vector<size_t> &cols,
                                   std::vector<double> &values ) const override;


    /** \brief  Given a row, retrieve the non-zero column indices of the matrix in compressed format
     * \param[in]  row Which row
     */
    std::vector<size_t> getColumnIDs( size_t row ) const override;

    virtual void getValuesByGlobalID( size_t num_rows,
                                      size_t num_cols,
                                      size_t *rows,
                                      size_t *cols,
                                      double *values ) const override;

    virtual void setScalar( double ) override;
    virtual void setDiagonal( Vector::const_shared_ptr in ) override;

    virtual void makeConsistent() override;
    virtual double L1Norm() const override;
    virtual Matrix::shared_ptr cloneMatrix() const override;
    virtual Vector::shared_ptr getRightVector() const override;
    virtual Vector::shared_ptr getLeftVector() const override;
    virtual Discretization::DOFManager::shared_ptr getRightDOFManager() const override;
    virtual Discretization::DOFManager::shared_ptr getLeftDOFManager() const override;
    virtual void fillComplete() override;
    virtual void setIdentity() override;
    virtual void zero() override;
};
} // namespace LinearAlgebra
} // namespace AMP

#include "ManagedEpetraMatrix.inline.h"
#endif
