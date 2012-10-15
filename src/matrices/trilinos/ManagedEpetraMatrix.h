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
      /** \brief  Parameters used to construct the matrix
        */
      boost::shared_ptr<ManagedEpetraMatrixParameters>         d_pParameters;
      
      /** \brief  Unimplemented constructor
        */
      ManagedEpetraMatrix();

      /** \brief  Unimplemented constructor
        */
      ManagedEpetraMatrix ( const ManagedEpetraMatrix &rhs );

      /** \brief  \f$A_{i,j}\f$ storage of off-core data
        */
      std::map<int,std::map<int,double> >  d_OtherData;

      /** \brief  Update data off-core
        */
      void  setOtherData ();

      virtual void multiply ( shared_ptr other_op , shared_ptr &result );

public:
      /** \brief Constructor
        * \param[in] p  The description of the matrix
        */
      ManagedEpetraMatrix( MatrixParameters::shared_ptr p );

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
      virtual Vector::shared_ptr  extractDiagonal ( Vector::shared_ptr buf = Vector::shared_ptr() );
      virtual void  scale ( double alpha );
      virtual void  axpy ( double alpha , const Matrix &rhs );
      virtual size_t  numGlobalRows () { return d_epetraMatrix->NumGlobalRows(); }
      virtual size_t  numGlobalColumns () { return d_epetraMatrix->NumGlobalCols(); }
      virtual void  addValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );
      virtual void  setValuesByGlobalID ( int   num_rows ,
                                          int   num_cols ,
                                          int  *rows ,
                                          int  *cols ,
                                          double  *values );

      virtual void  getRowByGlobalID ( int row ,
                                       std::vector<unsigned int> &cols,
                                       std::vector<double>       &values ) const;

      virtual void setScalar ( double );
      virtual void setDiagonal ( const Vector::shared_ptr &in );

      virtual void makeConsistent ();
      virtual double  L1Norm() const;
      virtual Matrix::shared_ptr cloneMatrix () const;
      virtual Vector::shared_ptr  getRightVector ();
      virtual Vector::shared_ptr  getLeftVector ();
      virtual Discretization::DOFManager::shared_ptr  getRightDOFManager ();
      virtual Discretization::DOFManager::shared_ptr  getLeftDOFManager ();
      virtual void fillComplete();
};


}
}

#include "ManagedEpetraMatrix.inline.h"
#endif


