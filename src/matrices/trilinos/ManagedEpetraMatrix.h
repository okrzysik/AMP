#ifndef included_AMP_ManagedEpetraMatrix
#define included_AMP_ManagedEpetraMatrix

#include <set>

#include "matrices/ManagedMatrix.h"
#include "matrices/trilinos/ManagedEpetraMatrixParameters.h"
#include "matrices/trilinos/EpetraMatrix.h"
#include "vectors/trilinos/EpetraVector.h"
#include "discretization/DOF_Manager.h"

#include <Epetra_FECrsMatrix.h>

namespace AMP {
namespace LinearAlgebra {


/** \class ManagedEpetraMatrix
  * \brief  A class that wraps an Epetra_CrsMatrix
  * \details  This class stores an Epetra_FECrsMatrix and provides
  * the AMP interface to this matrix.
  */
class ManagedEpetraMatrix : public EpetraMatrix ,
                              public ManagedMatrix
{
protected:
    //!  Parameters used to construct the matrix
    boost::shared_ptr<ManagedEpetraMatrixParameters>         d_pParameters;
      
    //!  Unimplemented constructor
    ManagedEpetraMatrix();

    //!  Unimplemented constructor
    ManagedEpetraMatrix ( const ManagedEpetraMatrix &rhs );

    //!  \f$A_{i,j}\f$ storage of off-core data
    std::map<int,std::map<int,double> >  d_OtherData;

    //!  Update data off-core
    void  setOtherData ();

    virtual void multiply ( shared_ptr other_op , shared_ptr &result );

public:
    /** \brief Constructor
      * \param[in] p  The description of the matrix
      */
    ManagedEpetraMatrix( boost::shared_ptr<ManagedEpetraMatrixParameters> p );

    /** \brief Constructor from Epetra_CrsMatrix
      * \param[in]  m  Matrix to wrap
      * \param[in]  dele  If true, this class deletes the matrix
      */
    ManagedEpetraMatrix ( Epetra_CrsMatrix *m , bool dele = false );

    /** \brief Destructor
      */
    virtual ~ManagedEpetraMatrix() {}

    virtual void  createValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );



    virtual void mult(const Vector::const_shared_ptr in, Vector::shared_ptr out);
    virtual void multTranspose (const Vector::const_shared_ptr in, Vector::shared_ptr out);
    virtual Vector::shared_ptr  extractDiagonal ( Vector::shared_ptr buf = Vector::shared_ptr() ) const;
    virtual void  scale ( double alpha );
    virtual void  axpy ( double alpha , const Matrix &rhs );
    virtual size_t  numGlobalRows () const { return d_epetraMatrix->NumGlobalRows(); }
    virtual size_t  numGlobalColumns () const { return d_epetraMatrix->NumGlobalCols(); }
    virtual void  addValuesByGlobalID ( int num_rows, int num_cols, int *rows, int *cols, double *values );
    virtual void  setValuesByGlobalID ( int num_rows, int num_cols, int *rows, int *cols, double *values );
    virtual void  getRowByGlobalID ( int row, std::vector<unsigned int> &cols, std::vector<double> &values ) const;
    virtual void  getValuesByGlobalID ( int num_rows, int num_cols, int *rows, int *cols, double *values ) const;

    virtual void setScalar ( double );
    virtual void setDiagonal ( Vector::const_shared_ptr in );

    virtual void makeConsistent ();
    virtual double  L1Norm() const;
    virtual Matrix::shared_ptr cloneMatrix () const;
    virtual Vector::shared_ptr  getRightVector () const;
    virtual Vector::shared_ptr  getLeftVector () const;
    virtual Discretization::DOFManager::shared_ptr  getRightDOFManager () const;
    virtual Discretization::DOFManager::shared_ptr  getLeftDOFManager () const;
    virtual void fillComplete();
    virtual void setIdentity();
    virtual void zero();
};


}
}

#include "ManagedEpetraMatrix.inline.h"
#endif


