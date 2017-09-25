#ifndef included_AMP_EpetraMatrix
#define included_AMP_EpetraMatrix

#include "matrices/Matrix.h"

#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {

/** \brief  Parameters class for creating an EpetraMatrix
 */
typedef MatrixParameters EpetraMatrixParameters;

/** \class EpetraMatrix
  * \brief A Matrix with an Epetra_CrsMatrix interface
  * \details  An EpetraMatrix presents an Epetra_Matrix class.
  * Given an AMP::LinearAlgebra::Matrix, this class can create an Epetra view
  * without copying the data.  As such, this class serves three
  * purposes:
  *  -# Provides an Epetra_CrsMatrix for derived classes to use, fill, manage, etc.
  *  -# Provides an interface for accessing this Epetra_CrsMatrix independent of base or derived
  classes
  *  -# Provides a static method for creating an Epetra_CrsMatrix view of an AMP vector.
  */

class EpetraMatrix : virtual public Matrix
{
private:
    EpetraMatrix();
    EpetraMatrix( const EpetraMatrix &rhs );
    explicit EpetraMatrix( MatrixParameters::shared_ptr );

protected:
    /** \brief Bare pointer to an Epetra_CrsMatrix
     */
    Epetra_CrsMatrix *d_epetraMatrix;

    /** \brief Range map for the Epetra_CrsMatrix
     */
    AMP::shared_ptr<Epetra_Map> d_RangeMap;

    /** \brief Domain map for the Epetra_CrsMatrix
     */
    AMP::shared_ptr<Epetra_Map> d_DomainMap;

    /** \brief Indicates if the destructor calls delete
     */
    bool d_DeleteMatrix;


    /** \brief Ensure Epetra methods return correctly
     * \param[in] err  The return value from the method
     * \param[in] func  The name of the Epetra method called
     * \details  Throws an execption if err != 0
     */
    void VerifyEpetraReturn( int err, const char *func ) const;

    /** \brief Constrcutor
     * \param[in]  m1  Rowmap to create the Epetra matrix
     * \param m2  Unused
     * \param[in] entriesRow  The number of entries in the matrix per local row
     */
    EpetraMatrix( Epetra_Map &m1, Epetra_Map *m2, int *entriesRow );

    /** \brief Constructor
     * \param[in] inMatrix  Matrix to wrap
     * \param[in] dele  If true, then this class will delete the Epetra_CrsMatrix
     */
    EpetraMatrix( Epetra_CrsMatrix *inMatrix, bool dele = false );


public:
    /** \brief Change the EpetraMaps for the matrix
     * \param[in] range  A vector that represents the range: y in y = A*x (row map)
     * \param[in] domain  A vector that represents the domain: x in y = A*x (column map)
     * \details  This does not change the matrix, just the maps stored above
     *
     */
    void setEpetraMaps( Vector::shared_ptr range, Vector::shared_ptr domain );

    /** \brief Destructor
     */
    virtual ~EpetraMatrix();

    /** \brief  Return an Epetra_CrsMatrix
     * \return An Epetra_CrsMatrix view of this matrix
     */
    virtual Epetra_CrsMatrix &getEpetra_CrsMatrix();

    /** \brief  Return an Epetra_CrsMatrix
     * \return An Epetra_CrsMatrix view of this matrix
     */
    virtual const Epetra_CrsMatrix &getEpetra_CrsMatrix() const;

    /** \brief  Create an EpetraMatrix view of an AMP::LinearAlgebra::Matrix
     * \param[in] p  The matrix to view
     * \return  An AMP:Matrix capable of casting to EpetraMatrix
     */
    static AMP::shared_ptr<EpetraMatrix> createView( shared_ptr p );

    /** \brief  A call-through to Epetra_CrsMatrix fillComplete
     */
    virtual void fillComplete();

    virtual Matrix::shared_ptr transpose() const;
};
} // namespace LinearAlgebra
} // namespace AMP

#include "EpetraMatrix.inline.h"

#endif
